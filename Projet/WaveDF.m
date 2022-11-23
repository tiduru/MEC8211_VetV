%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Projet
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cette fonction résout l'équation différentielle représentant l'équation 
% d'onde pour une corde tendue dans les conditions suivantes:
%
%   - Une corde de longueur L
%   - Un déplacement transversal initial nul aux frontières correspondant
%     à des condition frontière de Dirichlet, i.e. u(x=0,0)=u(x=L,0)=0
%   - Une vitesse initiale nulle, i.e. du/dt = 0 à t=0
%   - Un déplacement initial variable spécifié par une fonction f(x)
%
%  
% Variables
% ---------
%   entrée : Ntot   - Nombre de noeuds, Entier >= 3
%            dt     - Pas de temps [s], > 0
%            Ndt    - Nombre de pas de temps, Entier >= 1
%            L      - Longueur de la corde [m]
%            rho    - Densité linéaire de la corde [kg/m]
%            T      - Tension dans la corde [N]
%            b      - Damping coefficient [s^-1]
%           @fx     - Fonction retournant le déplacement initial [m]
%           @Sxt    - Fonction retournant le terme source [m/s^2]
%
%   sortie : u      - Déplacement transversal [m]. Taille Ndt+1 x Ntot
%                        Rangées (temps)  : Ndt + 1
%                        Colonnes (noeuds): Ntot
%                        Ex: u(1,1) = Déplacement initial au noeud 1.
%                            u(101,3) = Déplacement au noeud 3, temps 100*dt.
%
%            temps  - Temps discrets de la simulation [s]. Taille Ndt+1
%
% Historique
% 22-Nov-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, temps] = WaveDF(Ntot, dt, Ndt, L, rho, T, b, fx, Sxt)

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
if(L<=0)
   disp("La longueur de la corde tendue doit être plus grande que 0");
   valid = false;
end
if(rho <= 0)
   disp("La densité linéaire doit être pplus grande que 0 kg/m");
   valid = false;
end
if(T <= 0)
   disp("La tension dans la corde doit être plus grande que 0 N");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Intervalles dx
dx = L/(Ntot-1);

% Constantes
A = (1+b*dt/2)^-1;
B = (b*dt/2 - 1);
C = sqrt(T/rho)*dt/dx; % Nombre de Courant
if(C > 1)
   msg = sprintf("Le nombre de Courant de %.3f est supérieur à 1.", C);
   disp(msg);
   disp("Le schéma explicite peut être instable.");
   disp("Il est suggéré de diminuer dt ou Ntot.");
end
C2 = C^2; 

% Initialisation
temps = 0:dt:Ndt*dt;     % Temps discrets de la simulation [s]
u = zeros(Ndt, Ntot);    % Déplacements.

% Conditions frontière de Dirichlet
u(1,1) = 0;
u(1,Ntot) = 0;

% Déplacement initial à chaque noeud interne
f = zeros(Ntot);         
for i=2:Ntot-1
   f(i) = fx(dx*(i-1));
   u(1,i) = f(i);
end

% Équation spéciale pour le premier pas de temps
for i=2:Ntot-1
   Si = Sxt(dx*(i-1),0); % Terme source au noeud i, temps t=0
   u(2,i) = ((1-A*B)^-1)*(A*C2*(f(i+1)+f(i-1)) + 2*A*(1-C2)*f(i)+ A*Si*dt^2);
end

% Solution par différences finies (explicite)
for t=3:Ndt+1
  for i=2:Ntot-1
     u(t,i) = A*B*u(t-2,i) + ...
              A*C2*u(t-1,i+1) + ...
              2*A*(1-C2)*u(t-1,i) + ...
              A*C2*u(t-1,i-1) + ...
              A*(dt^2)*Sxt(dx*(i-1),dt*t);
  end
end
