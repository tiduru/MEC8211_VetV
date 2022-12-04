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
% est effectuée avec une solution manufacturée.
%
% Méthodologie: Le raffinement h est égal à dx. Puisque dt est lié à dx par
%               le nombre de Courant C, on choisit un C fixe de 0.9. On 
%               raffine ensuite h (divise par 2) à chaque itération et on
%               calcule le dt correspondant.
%
%
% Variables
% ---------
%   entrée : Ninit - Nombre de noeuds pour le maillage initial, Entier >= 3
%            Nraf  - Nombre de raffinements où dx est divisé par 2, Entier >= 1
%            tend  - Temps final de la simulation [s], > 0
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
% Exemple: 3 noeuds initial, 4 raffinements, temps final de 2 secondes
%          >> WaveVerifMMS(3, 4, 2)
%
% Historique
% 01-Dec-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WaveVerifMMS(Ninit, Nraf, tend)

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
if(tend <= 0)
   disp("Le temps final de l'analyse doit être supérieur à 0 secondes");
   valid = false;
end
if(~valid)
   error("Au moins un des arguments est invalide");
end

% Paramètres du problème
b = 0.5;           % Damping coefficient [s^-1]
L = 2.0;           % Longueur de la corde [m]
T = 39.66;         % Tension dans la corde [N]
rho = 50.33*10^-5; % Densité linéaire [kg/m]
c = sqrt(T/rho);   % Vitesse de propagation de l'onde
C = 0.9;           % Nombre de courant fixé à 0.9
k = 1;             % Paramètre variable
c2= T/rho;

% Polynome spatial et aussi conditions initiales 
fx = @(x) k*x^3-k*L*x^2;

% Solution manufacturée
u_ana = @(x,t) fx(x)*cos(c*t);

% Terme source obtenu de la solution manufacturée
Sxt = @(x,t) -b*c*sin(c*t)*fx(x)-c2*cos(c*t)*(fx(x) - 2*k*L+6*k*x);

% Erreurs L1, L2 et Linf
Niter = 1 + Nraf;
h = zeros(Niter, 1);
L1 = zeros(Niter, 1);
L2 = zeros(Niter, 1);
Linf = zeros(Niter, 2);

% Nombre de noeuds pour chaque itération
Ntot = Ninit;
for i=1:Niter 
   if(i > 1)
      Ntot = 1 + 2*(Ntot-1); % On double le nombre d'intervalles
   end
   h(i) = L/(Ntot-1);        % Intervalle pour cette itération [m]

   dt = (C/c)*h(i);
   Ndt = ceil(tend/dt);

   % Déplacements calculés par différences finies [m]
   [uO, t] = WaveDF(Ntot, dt, Ndt, L, rho, T, b, fx, Sxt); 
   udf = uO(Ndt+1,:);
   udf(end,:);

   % Déplacements analytiques [m]
   for j=1:Ntot
      uana(j) = u_ana((j-1)*h(i), Ndt*dt); 
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
   dh = log(h(p)/h(p+1));
   p_hat_L1(p)   = log(L1(p)/L1(p+1))/dh;
   p_hat_L2(p)   = log(L2(p)/L2(p+1))/dh;
   p_hat_Linf(p) = log(Linf(p)/Linf(p+1))/dh;
   disp(sprintf("pentes %d-%d: L1=%.2f, L2=%.2f, Linf=%.2f", p, p+1,...
        p_hat_L1(p), p_hat_L2(p), p_hat_Linf(p)));
end

% Création du graphe
figure
p1 = loglog(h, L1, '-s');
title(sprintf("Convergence Observée\n " + ...
   "Temps Final=%G s, Nombre de Courant C=%.1f", ...
   tend, C));
xlabel('h=\Delta{x}=\Delta{t}(v/C) [m]');
ylabel('Erreur [m]');
xlim padded
ylim padded
grid on
hold on
p2 = loglog(h, L2, '-s');
p3 = loglog(h, Linf, '-s');
p1(1).MarkerFaceColor = p1(1).Color;
p2(1).MarkerFaceColor = p2(1).Color;
p3(1).MarkerFaceColor = p3(1).Color;

% Triangle avec les pentes des points a et b (a > b)
a = 6; % 1er point
b = 7; % 2e point
triang_x = [h(b), h(a)];
triang_y = exp(interp1(log(h), log(L1), log(triang_x)));
triang_y = triang_y - 5e-4*[h(b)^p_hat_L1(a) h(a)^p_hat_L1(a)];
loglog(triang_x([1,2,2,1]), triang_y([1,1,2,1]), 'k');
text(h(a)+0.001, exp(0.5*(log(triang_y(1))+log(triang_y(2)))), ...
     sprintf('L_{1}=%.2f\nL_{2}=%.2f\nL_{\\infty}=%.2f', ...
     p_hat_L1(a), p_hat_L2(a),p_hat_Linf(a)))

end