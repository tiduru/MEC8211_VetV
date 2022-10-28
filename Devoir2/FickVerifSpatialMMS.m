%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 2
%
% Eduards Blandin     1893699
% Jacques Desfossés     61902
% Timothée Duruisseau 1949883
%
% Cet script permet la vérification de la convergence spatialle, par le
% biais de la fonction FickDFMMS, du schéma de différences finies
% implémenté dans la fonction FrickDF. La vérification sera effectuée par
% rapport à la solution manufacturée.
%
% Variables
% ---------
%   entrée : Ninit - Nombre de noeuds pour le maillage initial, Entier >= 3
%            Nraf  - Nombre de raffinements où h est divisé par 2, Entier >= 1
%            dt    - Pas de temps [an], > 0
%            Ndt   - Nombre de pas de temps, Entier >= 1
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
%   Exemple: 3 noeuds initial, 10 raffinements, 2 pas de temps de 0.1 an
%            directe: FickVerifSpatialMMS(3, 10, 0.1, 2)
%
% Historique
% 19-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FickVerifSpatialMMS(Ninit, Nraf, dt, Ndt)

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
t0 = 0.001; % Constante temporelle [an]
R = 0.5;    % Rayon du pilier [m]

% Erreurs L1, L2 et Linf
Niter = 1 + Nraf;
h  = zeros(Niter, 1);
L1 = zeros(Niter, 1);
L2 = zeros(Niter, 1);
Linf = zeros(Niter, 2);

 % Nombre de noeuds pour chaque itération
Ntot = Ninit;
for i=1:Niter 
   if(i > 1)
      Ntot = 1 + 2*(Ntot-1); % On double le nombre d'intervalles
   end
   h(i) = R/(Ntot-1);        % Intervalle pour cette itération [m]  

   % Concentrations calculées par différences finies [mol/m^3]
   [Cdf, Cmms] = FickDFMMS(Ntot, dt, Ndt, t0); 

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
p_hat_L1 = zeros(Niter-1);
p_hat_L2 = zeros(Niter-1);
p_hat_Linf = zeros(Niter-1);
for p=1:Niter-1
   dh = log(h(p)/h(p+1));
   p_hat_L1(p)   = log(L1(p)/L1(p+1))/dh;
   p_hat_L2(p)   = log(L2(p)/L2(p+1))/dh;
   p_hat_Linf(p) = log(Linf(p)/Linf(p+1))/dh;
   disp(sprintf("pentes %d-%d: L1=%f, L2=%f, Linf=%f", p, p+1,...
        p_hat_L1(p), p_hat_L2(p), p_hat_Linf(p)));
end

% Création du graphe
figure
p1 = loglog(h, L1, '-s');
title(sprintf("Convergence spatiale\n Ndt=%d, dt=%G ans", Ndt, dt));
xlabel('h [m]');
ylabel('Erreur [mol/m^{3}]');
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
a = 2; % 1er point
b = 3; % 2e point
triang_x = [h(b), h(a)];
triang_y = exp(interp1(log(h), log(L1), log(triang_x)));
triang_y = triang_y - [h(b)^p_hat_L1(a) h(a)^p_hat_L1(a)] ;
loglog(triang_x([1,2,2,1]), triang_y([1,1,2,1]), 'k');
text(h(a)+0.001, exp(0.5*(log(triang_y(1))+log(triang_y(2)))), ...
     sprintf('L_{1}=%.2f\nL_{2}=%.2f\nL_{\\infty}=%.2f', ...
     p_hat_L1(a), p_hat_L2(a),p_hat_Linf(a)))

% Légende
lgd = legend('L_{l}', 'L_{2}', 'L_{\infty}');
lgd.Location = 'northwest';

end