%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 1
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cet script permet la vérification des fonctions FrickDF et FrickDFStat
% qui calculent l'évolution de la concentration de sel à l'intérieur d'un
% pilier de béton à l'aide des différences finies. La vérification sera
% effectuée par rapport à la solution analytique stationnaire lorsque le
% terme source est constant.
%
% Variables
% ---------
%   entrée : Nmin      - Nombre de noeuds minimal, Entier >= 3
%            Nmax      - Nombre de noeuds maximal, Entier > Nmin
%            schema    - Schéma de différenciation: 1 - Ordre 1
%                                                   2 - Ordre 2
%            DFsol - Solution stationnaire par différences finies:
%                    'transitoire' - La solution transitoire roule
%                                    jusqu'à la solution stationnaire
%                    'directe'     - Résolution directe de l'équation
%                                    elliptique de la solution stationnaire
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
%   Exemple: Noeuds entre 3 et 10, schéma d'ordre 2, solution stationnaire
%            directe: FickVerifStat(3, 10, 2, "directe")
%
% Historique
% 13-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FickVerifStat(Nmin, Nmax, schema, DFsol)

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
if(schema ~=1 && schema ~=2)
   disp("L'ordre du schéma de différenciation doit être 1 ou 2");
   valid = false;
end
if(DFsol ~= "transitoire" && DFsol ~= "directe")
   disp("La solution stationnaire doit être transitoire ou directe");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Paramètres du problème
R = 0.5;  % Rayon du pilier [m]

% Paramètres pour la méthode "transitoire"
dt     = 1;    % Pas de temps [an]
Ndt    = 1000; % Nombre de pas de temps
tsMeth = 0;    % Méthode "constante" pour le terme source

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
   Cdf = zeros(Ntot, 1);
   if (DFsol == "transitoire")
      [CO, t] = FickDF(Ntot, dt, Ndt, schema, tsMeth); 
      Cdf = CO(Ndt+1,:);
   elseif (DFsol == "directe")
      Cdf = FickDFStat(Ntot, schema);
   end

   % Concentrations analytiques [mol/m^3]
   Cana = FickAnaStat(Ntot); 

   % Calcul des erreurs
   for j=1:Ntot
      diffabs = abs(Cdf(j) - Cana(j));
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
disp(sprintf("pentes O(%d): L1=%f, L2=%f, Linf=%f", ...
     schema, p_hat_L1, p_hat_L2, p_hat_Linf))

% Création du graphe
figure
loglog(h, L1, h, L2, h, Linf);
title(sprintf("Erreurs pour un schéma de différenciation O(%d)\n Méthode %s", ...
              schema, DFsol));
xlabel('h [m]');
ylabel('Erreur');
legend('L_{l}', 'L_{2}', 'L_{\infty}');
grid on

end
