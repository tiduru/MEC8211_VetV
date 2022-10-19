%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 2
%
% Eduards Blandin     1893699
% Jacques Desfossés     61902
% Timothée Duruisseau 1949883
%
% Cette fonction résout l'équation différentielle représentant la 2e loi
% de Frick exprimée en coordonnées cylindriques avec un terme source
% associé à la solution manufacturée suivante:
%
% C_hat(t,r)=(1-1/(1+t))*(cos(pi*r^2/(2*R^2)))^2
% 
% La fonction calcule C(r,t) pour:
%
%   - Un pilier de béton de rayon R submergé dans une eau saline
%   - Un pilier (cylindre) infiniment haut
%   - Une concentration nulle à la surface du pilier correspondant à une
%     condition frontière de Dirichlet, i.e. C(t,r=R) = 0
%   - Un flux nul en r=0 correspondant à une condition frontière de Neumann
%   - Une concentration initiale nulle dans le pilier C(t=0, r) = 0
%
%  
% Variables
% ---------
%   entrée : Ntot   - Nombre de noeuds, Entier >= 3
%            dt     - Pas de temps [an], > 0
%            Ndt    - Nombre de pas de temps, Entier >= 1
%
%   sortie : C      - Concentrations [mol/m^3]. Taille Ntot
%            Cmms   - Concentrations de la solution manufacturée
%                     [mol/m^3]. Taille Ntot
%            rpos   - Position des noeuds (r_i). Taille Ntot
%
%   test : 50 noeuds, 30 incréments d'un mois:
%          [C, Cmms, r] = FrickDFMMS(50, 1/12, 30);
%
% Historique
% 19-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, Cmms, rpos] = FickDFMMS(Ntot, dt, Ndt)

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
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Données du problème
R    = 0.5;    % Rayon du pilier de béton [m]
Deff = 1E-10;  % Cefficient de diffusion effectif du sel [mol/m^3/s]
k    = 4E-9;   % Constante de réaction [1/s]

% On exprime Deff, k et Scst en [1/an]
secParAn = 31536000;
Deff = Deff * secParAn; % [mol/m^3/an]
k = k * secParAn;       % [1/an]

% Intervalles h
h = R/(Ntot-1);

% Initialisation
temps = 0:dt:Ndt*dt;     % Temps discrets de la simulation [an]
rpos = 0:h:R;            % Position radiale des noeuds [m]
C = zeros(Ndt, Ntot);    % Concentrations. C(t=0,:)=0 -> Condition init.
Cana = zeros(Ndt, Ntot); % Concentrations exactes pour la solution manuf.
M = zeros(Ntot, Ntot);   % Matrice [M] à gauche du système linéaire
V = zeros(Ntot, 1);      % Vecteur {V} à droite du système linéaire

% Terme source S(t,r)
syms t;
syms r;
syms Rs;
syms ks;
syms Ds;
syms f(t);
syms A(r,Rs);
syms C_hat(t,r,Rs);          % Solution manufacturée
syms dC_hat_sur_rdr(t,r,Rs); % (1/r) * d(C_hat)/dr
syms S(t,r,Rs,ks,Ds);        % Terme source

f(t) = 1-1/(1+t);
A(r,Rs) = (pi/2)*(r/Rs)^2; 
C_hat(t,r,Rs)=f(t)*cos(A(r,Rs))^2; % Solution manufacturée
dC_hat_sur_rdr(t,r,Rs) = -f(t)*(pi/Rs^2)*sin(2*A(r,Rs));
S(t,r,Rs,ks,Ds) = diff(C_hat(t,r,Rs),t,1) - Ds*dC_hat_sur_rdr(t,r,Rs) - ...
                  Ds*diff(C_hat(t,r,Rs),r,2) + ks*C_hat(t,r,Rs);

% Valeurs de [M] et {V} constantes pour tous les pas de temps dû aux 
% conditions frontières

% Neumann - Noeud à r=0, flux nul - de Gear avant (ordre 2)
M(1,1:3) = [-3 4 -1];
V(1,1) = 0;

% Dirichlet - Noeud à r=R, concentration constante C(t,r=R)=0
M(Ntot,Ntot) = 1;
V(Ntot,1)    = 0;

% Valeur de E constante pour toute la simulation
E = Deff*h*dt;

% Solution par différences finies
for t=2:Ndt+1
   % Système matriciel [M]{C} = {V}
   for i=2:Ntot-1
      Si = S(temps(t),rpos(i),R,k,Deff);
      Ai = rpos(i)*Deff*dt;
      Bi = rpos(i)*(h^2);
      Fi = Bi*C(t-1,i) + Si*Bi*dt;
      ki = k*Bi*dt;

      M(i,i-1) = 0.5*E-Ai;
      M(i,i)   = Bi + 2*Ai + ki;
      M(i,i+1) = -0.5*E - Ai;

      V(i,1) = Fi;
   end

   % Solution au temps t: [C]=[M]^-1 {V}. On laisse Matlab déterminer le
   % meilleur algorithme en utilisant l'opérateur '\'
   C(t,:) = M\V;
end

% On ne conserve que les concentrations du dernier pas de temps
C = C(end,:)';

% Concentrations exactes de la solution manufacturées
Cmms = zeros(Ntot,1);
for i=1:Ntot
   Cmms(i) = C_hat(temps(end),rpos(i),R);
end
