%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 1
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cette fonction résout l'équation différentielle elliptique représentant
% la 2e loi de Frick exprimée en coordonnées cylindriques lorsque le régime
% stationnaire est atteint pour un terme source constant. La concentration
% C(r) est calculée pour:
%
%   - Un pilier de béton de rayon R submergé dans une eau saline
%   - Un pilier (cylindre) infiniment haut
%   - Une concentration Ce=10 mol/m^3 constante à la surface du pilier
%     correspondant à une condition frontière de Dirichlet
%   - Un flux nul en r=0 correspondant à une condition frontière de Neumann
%   - Une concentration initiale nulle dans le pilier C(r<R,0) = 0
%
%  
% Variables
% ---------
%   entrée : Ntot   - Nombre de noeuds, Entier >= 3
%            schema - Schéma de différenciation: 1 - Ordre 1
%                                                2 - Ordre 2
%
%   sortie : C      - Concentrations [mol/m^3]
%                     rangées (noeuds): Ntot
%                     Ex: C(1) = Concentration au noeud 1.
%                         C(3) = Concentration au noeud 3.
%
%   test : 50 noeuds, schéma d'ordre 2: C = FrickDF(50, 2);
%
% Historique
% 11-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = FickDFStat(Ntot, schema)

% Validation de toutes les variables d'entrée
valid = true;
if(Ntot < 3 || (floor(Ntot) ~= ceil(Ntot)))
   disp("Le nombre de noeuds Ntot doit être un entier >= 3");
   valid = false;
end
if(schema ~=1 && schema ~=2)
   disp("L'ordre du schéma de différenciation doit être 1 ou 2");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Données du problème
R    = 0.5;    % Rayon du pilier de béton [m]
Ce   = 10;     % Concentration à la surface du pilier [mol/m^3]
Deff = 10E-10; % Cefficient de diffusion effectif du sel [mol/m^3/s]
Scst = 1E-8;   % Terme source constant pour solution analytique [mol/m^3/s]

% On exprime Deff et Scst en [1/an]
secParAn = 31536000;
Deff = Deff * secParAn; % [mol/m^3/an]
Scst = Scst * secParAn; % [mol/m^3/an]

% Intervalles h
h = R/(Ntot-1);

% Initialisation
C = zeros(Ntot, 1);     % Concentrations
M = zeros(Ntot, Ntot);  % Matrice [M] à gauche du système linéaire
V = zeros(Ntot, 1);     % Vecteur {V} à droite du système linéaire

% Valeurs de [M] et {V} constantes dû aux conditions frontières

% Neumann - Noeud à r=0, flux nul
if(schema == 1)
   M(1,1:2) = [1 -1]; % dérivée première avant (ordre 1)
elseif(schema == 2)
   M(1,1:3) = [-3 4 -1]; % de Gear avant (ordre 2)
end
V(1,1) = 0;

% Dirichlet - Noeud à r=R, concentration constante Ce
M(Ntot,Ntot) = 1;
V(Ntot,1)    = Ce;

% Système matriciel [M]{C} = {V}
for i=2:Ntot-1
   ri = (i-1)*h; % Position radiale du noeud 

   if(schema == 1)
      M(i,i-1) = ri;
      M(i,i)   = -(2*ri + h);
      M(i,i+1) = ri + h;
   elseif(schema == 2)
      M(i,i-1) = ri-0.5*h;
      M(i,i)   = -2*ri;
      M(i,i+1) = ri + 0.5*h;
   end
   V(i,1) = (Scst*(h^2)*ri)/Deff;
end

% Solution au temps t: [C]=[M]^-1 {V}. On laisse Matlab déterminer le
% meilleur algorithme en utilisant l'opérateur '\'
C = M\V;

end
