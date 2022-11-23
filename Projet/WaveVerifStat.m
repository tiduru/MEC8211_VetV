%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Projet
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cet script permet la vérification de la fonctions WaveDF qui calcule le
% déplacement transversal d'une onde dans une corde tendue. La vérification
% est effectuée par rapport à la solution analytique obtenue lorsque le
% déplacement initial est f(x) = sin(pi*x/L), l'amortissement b est nul et
% le terme source est nul.
%
% Variables
% ---------
%   entrée : Ninit - Nombre de noeuds pour le maillage initial, Entier >= 3
%            Nraf  - Nombre de raffinements où dx est divisé par 2, Entier >= 1
%            dt    - Pas de temps [s], > 0
%            Ndt   - Nombre de pas de temps, Entier >= 1
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
% Exemple: 8 noeuds initial, 4 raffinements, 10,000 pas de temps de 100 ms
%          >> WaveVerifStat(8, 4, 1e-5,100000)
%
% Historique
% 22-Nov-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WaveVerifStat(Ninit, Nraf, dt, Ndt)

% Validation de toutes les variables d'entrée
valid = true;
if(Ninit < 3 || (floor(Ninit) ~= ceil(Ninit)))
   disp("Le nombre de noeuds initial Ninit doit être un entier >= 3");
   valid = false;
end
if(Nraf < 1 || (floor(Nraf) ~= ceil(Nraf)))
   disp("Le nombre de raffinements Nraf doit être un entier >= 1");
   valid = false;
end
if(dt <= 0)
   disp("L'intervalle de temps dt doit être > 0");
   valid = false;
end
if(Ndt < 1 || (floor(Ndt) ~= ceil(Ndt)))
   disp("Le nombre de pas de temps Ndt doit être une entier >= 1");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Paramètres du problème
b = 0;             % Damping coefficient [s^-1]
L = 2.0;           % Longueur de la corde [m]
T = 39.66;         % Tension dans la corde [N]
rho = 50.33*10^-5; % Densité linéaire [kg/m]
c = sqrt(T/rho);   % Vitesse de propagation de l'onde

% Conditions initiales 
fx = @(x) sin(pi*x/L);

% Terme source
Sxt = @(x,t) 0;

% Solution analytique
u_ana = @(x,t) sin(pi*x/L)*cos(pi*c*t/L);

% Paramètres transitoires
%dt     = 0.00001;  % Pas de temps [s]
%Ndt    = 20000;    % Nombre de pas de temps

% Erreurs L1, L2 et Linf pour un nombre de noeuds dans le
Niter = 1 + Nraf;
dx    = zeros(Niter, 1);
L1 = zeros(Niter, 1);
L2 = zeros(Niter, 1);
Linf = zeros(Niter, 2);

% Nombre de noeuds pour chaque itération
Ntot = Ninit;
for i=1:Niter 
   if(i > 1)
      Ntot = 1 + 2*(Ntot-1); % On double le nombre d'intervalles
   end
   dx(i) = L/(Ntot-1);       % Intervalle pour cette itération [m]    

   % Déplacements calculés par différences finies [m]
   [uO, t] = WaveDF(Ntot, dt, Ndt, L, rho, T, b, fx, Sxt); 
   udf = uO(Ndt+1,:);
   udf(end,:);

   % Déplacements analytiques [m]
   for j=1:Ntot
      uana(j) = u_ana((j-1)*dx(i), Ndt*dt); 
   end

   % Calcul des erreurs
   for j=1:Ntot
      diffabs = abs(udf(j) - uana(j));
      L1(i) = L1(i) + diffabs;
      L2(i) = L2(i) + diffabs^2;
      if(diffabs > Linf(i))
         Linf(i) = diffabs;
      end
   end
   L1(i) = L1(i)/Ntot;
   L2(i) = sqrt(L2(i)/Ntot);  
end

% Ordre de précision observé
p_hat_L1 = zeros(Niter-1);
p_hat_L2 = zeros(Niter-1);
p_hat_Linf = zeros(Niter-1);
for p=1:Niter-1
   dh = log(dx(p)/dx(p+1));
   p_hat_L1(p)   = log(L1(p)/L1(p+1))/dh;
   p_hat_L2(p)   = log(L2(p)/L2(p+1))/dh;
   p_hat_Linf(p) = log(Linf(p)/Linf(p+1))/dh;
   disp(sprintf("pentes %d-%d: L1=%.2f, L2=%.2f, Linf=%.2f", p, p+1,...
        p_hat_L1(p), p_hat_L2(p), p_hat_Linf(p)));
end

% Création du graphe
figure
p1 = loglog(dx, L1, '-s');
title(sprintf("Convergence spatiale\n Ndt=%d, dt=%G s", Ndt, dt));
xlabel('dx [m]');
ylabel('Erreur [m]');
xlim padded
ylim padded
grid on
hold on
p2 = loglog(dx, L2, '-s');
p3 = loglog(dx, Linf, '-s');
p1(1).MarkerFaceColor = p1(1).Color;
p2(1).MarkerFaceColor = p2(1).Color;
p3(1).MarkerFaceColor = p3(1).Color;

% Triangle avec les pentes des points a et b (a > b)
a = Niter-1; % 1er point
b = Niter;   % 2e point
triang_x = [dx(b), dx(a)];
triang_y = exp(interp1(log(dx), log(L1), log(triang_x)));
triang_y = triang_y - 10*[dx(b)^p_hat_L1(a) dx(a)^p_hat_L1(a)];
loglog(triang_x([1,2,2,1]), triang_y([1,1,2,1]), 'k');
text(dx(a)+0.001, exp(0.5*(log(triang_y(1))+log(triang_y(2)))), ...
     sprintf('L_{1}=%.2f\nL_{2}=%.2f\nL_{\\infty}=%.2f', ...
     p_hat_L1(a), p_hat_L2(a),p_hat_Linf(a)))

% Légende
lgd = legend('L_{l}', 'L_{2}', 'L_{\infty}');
lgd.Location = 'northwest';

end