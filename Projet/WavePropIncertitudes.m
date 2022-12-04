%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Projet
%
% Eduards Blandin
% Jacques Desfossés
% Timothée Duruisseau
%
% Cet script effectue la propagation des incertitudes pour le problème
% de la propagation d'une onde sur une corde tendue pincée. 
%
% La propagation est effectuée pour une incertitude épistémique sur
% l'amortissement b et pour des incertitudes aléatoires gaussiennes pour
% la tension T et la densité linéaire rho.
%
%   sortie : 1) CDF des incertitudes
%            2) Matrice umn sauvegardée pour refaire les graphes au besoin
%
% Historique
% 02-Déc-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WavePropIncertitudes()

% Valeurs constantes du problème
L = 2.0;           % Longueur de la corde [m]
g = 9.81;          % Accélération gravitationnelle [m/s^2]

% Incertitudes
bl = 0.1;          % Valeur minimale du coefficient d'amortissement [s^-1]
bh = 0.5;          % Valeur maximale du coefficient d'amortissemnet [s^-1]

muT = 39.66;       % Valeur moyenne de la tension dans la corde [N]
stdT = 3.0;        % Écart-type de la tension dans la corde [N]

murho = 50.33E-5;  % Valeur moyenne de la densité linéaire [kg/m]
stdrho = 5E-5;     % Écart-type de la densité linéaire [kg/m]

% Conditions initiales (corde pincée en x=L/2, hauteur h)
h = 0.5;           % Hauteur pincée [m]
fx = @(x) (2*h*x/L)*(x<=L/2)+(2*h*(1-x/L))*(x>L/2);

% Terme source (gravité constante)
Sxt = @(x,t) -g;

% Paramètres de discrétisation/solution transitoire
Ntot = 101;          % Nombre de noeuds total
tend = 1;            % Temps final de l'analyse [s]
dt = 1E-5;           % Pas de temps [s]
Nc = (Ntot-1)/2;     % Noeud central
Ndt = ceil(tend/dt); % Nombre de pas de temps

% Paramètres Monte Carlo
M = 10;   % Nombres d'échantillons pour l'incertitude épistémique sur "b"
N = 3000; % Nombres d'échantillons pour l'incertitude aléatoire sur T et rho

% Échantillons d'amortissement
seedb = 259;
rng(seedb);
b = sort(unifrnd(bl, bh, [1,M])); % Courbes tracées en ordre croissant de b

% Échantillons de tension (distribution gaussienne)
seedT = 123;
rng(seedT);
T = normrnd(muT, stdT, [1,N]);

% Échantillons de densité linéaire (distribution gaussienne)
seedrho = 543; % seed pour la gération des échantilles de densité linéaire
rng(seedrho);
rho = normrnd(murho, stdrho, [1,N]);

% Simulation Monte Carlo - Calcul de la SRQ
umn = zeros(M,N);
for m=1:M
   m % La simulation est longue, on imprime le no. d'échantillon m
   for n=1:N
      % Déplacements calculés par différences finies [m]
      [udf, t] = WaveDF(Ntot, dt, Ndt, L, rho(n), T(n), b(m), fx, Sxt);
      
      % SRQ - Amplitude du déplacement au noeud central à t=tend=1 s
      umn(m,n) = abs(udf(end,Nc));
   end
end

% On sauvegarde les matrices pour refaire les graphes si nécessaire
save('umn.mat', 'umn')

% Superposition des CDF pour chaque valeur de "b"
lgd = cell(M,1);
figure
hold on

for m=1:M
   ecdf(umn(m,:));
   lgd{m} = sprintf('b=%.3f', b(m));
end

ul=min(umn,[],'all');
uh=max(umn,[],'all');
axis([ul-.1*abs(ul) uh+0.1*abs(uh) 0 1]);
lgdn = legend(lgd);
lgdn.Location = 'northwest';
title('Fonction de distribution cumulative (CDF)');
xlabel('|u_{c}| [m]')
ylabel('Probabilité cumulée')

end
