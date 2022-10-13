%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 1
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cette fonction résout l'équation différentielle représentant la 2e loi
% de Frick exprimée en coordonnées cylindriques et calcule C(r,t) pour:
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
%            dt     - Pas de temps [an], > 0
%            Ndt    - Nombre de pas de temps, Entier >= 1
%            schema - Schéma de différentiation: 1 - Ordre 1
%                                                2 - Ordre 2
%            tsMeth - Terme source: 0 - Constant S=1E-8 [mol/m^3/an]
%                                   1 - 1er ordre, S=kC avec k=4E-9 [an^-1]
%
%   sortie : C      - Concentrations [mol/m^3]. Taille Ndt+1 x Ntot
%                     Rangées (temps)  : Ndt + 1
%                     Colonnes (noeuds): Ntot
%                     Ex: C(1,1) = Concentration initiale au noeud 1.
%                         C(101,3) = Concentration au noeud 3, temps 100*dt.
%
%            temps  - Temps discrets de la simulation [an]. Taille Ndt+1
%
%   test : 50 noeuds, 30 incréments d'une année, schéma d'ordre 2,
%          terme src cst: C = FrickDF(50, 1, 30, 2, 0);
%
% Historique
% 02-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, temps] = FickDF(Ntot, dt, Ndt, schema, tsMeth)

% Validation de toutes les variables d'entrée
valid = true;
if(Ntot < 3 || (floor(Ntot) ~= ceil(Ntot)))
   disp("Le nombre de noeuds Ntot doit être un entier >= 3");
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
if(schema ~=1 && schema ~=2)
   disp("L'ordre du schéma de différenciation doit être 1 ou 2");
   valid = false;
end
if(tsMeth ~=0 && tsMeth ~=1)
   disp("La méthode pour le terme source doit être 0 (cst) ou 1 (1er ordre");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Données du problème
R    = 0.5;    % Rayon du pilier de béton [m]
Ce   = 10;     % Concentration à la surface du pilier [mol/m^3]
Deff = 10E-10; % Cefficient de diffusion effectif du sel [mol/m^3/s]
k    = 4E-9;   % Constante de réaction [1/s]
Scst = 1E-8;   % Terme source constant pour solution analytique [mol/m^3/s]

% On exprime Deff, k et Scst en [1/an]
secParAn = 31536000;
Deff = Deff * secParAn; % [mol/m^3/an]
k = k * secParAn;       % [1/an]
Scst = Scst * secParAn; % [mol/m^3/an]

% Intervalles h
h = R/(Ntot-1);

% Initialisation
temps = 0:dt:Ndt*dt;     % Temps discrets de la simulation [an]
C = zeros(Ndt, Ntot);    % Concentrations. C(t=0,:)=0 -> Condition init.
M = zeros(Ntot, Ntot);   % Matrice [M] à gauche du système linéaire
V = zeros(Ntot, 1);      % Vecteur {V} à droite du système linéaire

% Ajustement de S et k en fonction de la méthode pour le terme source
if (tsMeth == 0) % Terme source constant
   k = 0;
elseif (tsMeth == 1) % Terme source de 1er ordre
   Scst = 0;
end

% Valeurs de [M] et {V} constantes pour tous les pas de temps dû aux 
% conditions frontières

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

% Valeur de E constante pour toute la simulation
E = Deff*h*dt;

% Solution par différences finies
for t=2:Ndt+1
   % Système matriciel [M]{C} = {V}
   for i=2:Ntot-1
      ri = (i-1)*h; % Position radiale du noeud
      Ai = ri*Deff*dt;
      Bi = ri*(h^2);
      Fi = Bi*C(t-1,i) - Scst*Bi*dt;
      ki = k*Bi*dt;

      if(schema == 1)
         M(i,i-1) = -Ai;
         M(i,i)   = Bi + 2*Ai + E + ki;
         M(i,i+1) = -E - Ai;
      elseif(schema == 2)
         M(i,i-1) = 0.5*E-Ai;
         M(i,i)   = Bi + 2*Ai + ki;
         M(i,i+1) = -0.5*E - Ai;
      end
      V(i,1) = Fi;
   end

   % Solution au temps t: [C]=[M]^-1 {V}. On laisse Matlab déterminer le
   % meilleur algorithme en utilisant l'opérateur '\'
   C(t,:) = M\V;
end
