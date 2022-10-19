%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 2
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cet script permet la vérification, par le biais de la fonction FickDFMMS,
% du schéma de différences finies implémenté dans la fonction FrickDF.
% La vérification sera effectuée par rapport à la solution manufacturée.
%
% Variables
% ---------
%   entrée : Nmin   - Nombre de noeuds minimal, Entier >= 3
%            Nmax   - Nombre de noeuds maximal, Entier > Nmin
%            dt     - Pas de temps [an], > 0
%            Ndt    - Nombre de pas de temps, Entier >= 1
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
%   Exemple: Noeuds entre 3 et 10, schéma d'ordre 2, solution stationnaire
%            directe: FickVerifStat(3, 10, 2, "directe")
%
% Historique
% 19-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FickVerifMMS(Nmin, Nmax, dt, Ndt)

% Validation de toutes les variables d'entrée
valid = true;
if(Nmin < 3 || (floor(Nmin) ~= ceil(Nmin)))
   disp("Le nombre de noeuds minimal Nmin doit être un entier >= 3");
   valid = false;
end
if(Nmax <= Nmin || (floor(Nmax) ~= ceil(Nmax)))
   disp("Le nombre de noeuds maximal Nmax doit être un entier > Nmin");
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
R = 0.5;  % Rayon du pilier [m]

% Erreurs L1, L2 et Linf pour un nombre de noeuds dans le
% maillage: Nmin <= Ntot <= Nmax
Niter = Nmax - Nmin + 1;
h  = zeros(Niter, 1);
L1 = zeros(Niter, 1);
L2 = zeros(Niter, 1);
Linf = zeros(Niter, 2);

for i=1:Niter
   Ntot = Nmin + (i-1);  % Nombre de noeuds pour cette itération
   h(i) = R/(Ntot-1);    % Intervalle pour cette itération [m]  

   % Concentrations calculées par différences finies [mol/m^3]
   [Cdf, Cmms] = FickDFMMS(Ntot, dt, Ndt); 

   % Calcul des erreurs
   for j=1:Ntot
      diffabs = abs(Cdf(j) - Cmms(j));
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
a = 1; % 1er point pour déterminer la pente
b = 2; % 2e  point pour déterminer la pente
dh = log(h(a)/h(b));

p_hat_L1   = log(L1(a)/L1(b))/dh;
p_hat_L2   = log(L2(a)/L2(b))/dh;
p_hat_Linf = log(Linf(a)/Linf(b))/dh;
disp(sprintf("pentes: L1=%f, L2=%f, Linf=%f", p_hat_L1, p_hat_L2, p_hat_Linf));

% Création du graphe
figure
loglog(h, L1, h, L2, h, Linf);
title(sprintf("Erreurs après %d pas de temps de %G an", ...
              Ndt, dt));
xlabel('h [m]');
ylabel('Erreur [m]');
legend('L_{l}', 'L_{2}', 'L_{\infty}');
grid on

end
