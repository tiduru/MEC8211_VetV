%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 1
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cet script permet la vérification de la fonction FrickDF qui calcule 
% l'évolution de la concentration de sel à l'intérieur d'un pilier de
% béton à l'aide des différences finies. La vérification sera effectuée par
% rapport à la solution analytique stationnaire lorsque le terme source est
% constant.
%
% Historique
% 03-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

% Paramètres de la solution par différences finies
Nmin   = 3;    % Nombre de noeuds minimal
Nmax   = 10;   % Nombre de noeuds maximal
R      = 0.5;  % Rayon du pilier [m]
dt     = 1;    % Pas de temps [an]
Ndt    = 1000; % Nombre de pas de temps
tsMeth = 0;    % Méthode "constante" pour le terme source

% Calcul des erreurs L1, L2 et Linf pour un nombre de noeuds dans le
% maillage: Nmin <= Ntot <= Nmax
Niter = Nmax - Nmin + 1;
h  = zeros(Niter, 1);
L1 = zeros(Niter, 2);
L2 = zeros(Niter, 2);
Linf = zeros(Niter, 2);
for i=1:Niter
   Ntot = Nmin + (i-1);  % Nombre de noeuds pour cette itération
   h(i) = R/(Ntot-1);    % Intervalle pour cette itération [m]
   %[CO1, t1] = FickDF(Ntot, dt, Ndt, 1, tsMeth); % Concentrations, schéma O(1)
   %[CO2, t2] = FickDF(Ntot, dt, Ndt, 2, tsMeth); % Concentrations, schéma O(2)
   CO1 = FickDFStat(Ntot, 1);
   CO2 = FickDFStat(Ntot, 2);
   Cana = FickAnaStat(Ntot); % Concentrations analytiques [mol/m^3]

   for j=1:Ntot
      % Ordre 1
      %diffabs1 = abs(CO1(Ndt+1,j) - Cana(j));
      diffabs1 = abs(CO1(j) - Cana(j));
      L1(i,1) = L1(i,1) + diffabs1;
      L2(i,1) = L2(i,1) + diffabs1^2;
      if(diffabs1 > Linf(i,1))
         Linf(i,1) = diffabs1;
      end
      
      % Ordre 2
      %diffabs2 = abs(CO2(Ndt+1,j) - Cana(j));    
      diffabs2 = abs(CO2(j) - Cana(j));    
      L1(i,2) = L1(i,2) + diffabs2;
      L2(i,2) = L2(i,2) + diffabs2^2;
      if(diffabs2 > Linf(i,2))
         Linf(i,2) = diffabs2;
      end
   end

   % Ordre 1
   L1(i,1) = L1(i,1)/Ntot;
   L2(i,1) = sqrt(L2(i,1)/Ntot);
   
   % Ordre 2
   L1(i,2) = L1(i,2)/Ntot;
   L2(i,2) = sqrt(L2(i,2)/Ntot);
end

% Ordre de précision observé
a = 1; % 1er point pour déterminer la pente
b = 3; % 2e  point pour déterminer la pente
dh = log(h(a)/h(b));

p_hat_L1   = log(L1(a,1)/L1(b,1))/dh;
p_hat_L2   = log(L2(a,1)/L2(b,1))/dh;
p_hat_Linf = log(Linf(a,1)/Linf(b,1))/dh;
disp(sprintf("p_hat O(1): L1=%f, L2=%f, Linf=%f", p_hat_L1, p_hat_L2, p_hat_Linf))

p_hat_L1   = log(L1(a,2)/L1(b,2))/dh;
p_hat_L2   = log(L2(a,2)/L2(b,2))/dh;
p_hat_Linf = log(Linf(a,2)/Linf(b,2))/dh;
disp(sprintf("p_hat O(2): L1=%f, L2=%f, Linf=%f", p_hat_L1, p_hat_L2, p_hat_Linf))


% Création des graphes
figure
subplot(2,1,1);
loglog(h, L1(:,1), h, L2(:,1), h, Linf(:,1));
title('Erreurs pour un schéma de différenciation O(1)');
xlabel('h');
ylabel('Erreur');
legend('L_{l}', 'L_{2}', 'L_{\infty}');

subplot(2,1,2);
loglog(h, L1(:,2), h, L2(:,2), h, Linf(:,2));
title('Erreurs pour un schéma de différenciation O(2)');
xlabel('h');
ylabel('Erreur');
legend('L_{l}', 'L_{2}', 'L_{\infty}');