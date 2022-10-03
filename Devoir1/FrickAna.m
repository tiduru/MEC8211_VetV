%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 1
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cette fonction calcule la concentration analytique de sel dans le pilier
% de béton lorsque le terme source est constant et que le régime
% stationnaire est atteint.
%
%  
% Variables
% ---------
%   entrée : Ntot   - Nombre de noeuds, Entier >= 3
%
%   sortie : C      - Concentrations [mol/m^3]
%
% Historique
% 03-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cana = FrickAna(Ntot)

% Validation de la variable d'entrée
if(Ntot < 3 || (floor(Ntot) ~= ceil(Ntot)))
   error("Le nombre de noeuds Ntot doit être un entier >= 3");
   valid = false;
end

% Données du problème
R    = 1;      % Rayon du pilier de béton [m]
Ce   = 10;     % Concentration à la surface du pilier [mol/m^3]
Deff = 10E-10; % Cefficient de diffusion effectif du sel [mol/m^3/an]
Scst = 1E-8;   % Terme source constant pour solution analytique [mol/m^3/an]

% Correction temporaire assumant que les données Deff, k et Scst sont en
% 1/s en non en 1/an. On les convertit donc en 1/an.
secParAn = 365*24*60*60;
Deff = Deff * secParAn;
Scst = Scst * secParAn;

% Intervalles h
h = R/(Ntot-1);

% Calcul de la concentration analytique
Cana = zeros(Ntot,1);
for i=1:Ntot
   ri = (i-1)*h;
   Cana(i) = 0.25*(Scst/Deff)*(ri^2 - R^2) + Ce;
end
