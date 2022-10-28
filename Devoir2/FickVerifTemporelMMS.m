%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% MEC8211 A2022 Devoir 2
%
% Eduards Blandin     1893699
% Jacques Desfossés     61902
% Timothée Duruisseau 1949883
%
% Cet script permet la vérification de la convergence temporelle, par le
% biais de la fonction FickDFMMS, du schéma de différences finies
% implémenté dans la fonction FickDF. La vérification sera effectuée
% par rapport à la solution manufacturée.
%
% Variables
% ---------
%   entrée : Ntot    - Nombre de noeuds pour le maillage, Entier >= 3
%            Nraf    - Nombre de raffinements où dt est divisé par 2,
%                      Entier >= 1
%            dtInit  - Pas de temps initial [an], > 0
%            NdtInit - Nombre de pas de temps initial, Entier >= 1
%
%   sortie : 1) Graphe des erreurs L1, L2 et Linfini
%            2) Impression des pentes L1, L2 et Linfini
%
%   Exemple: 10 noeuds, 10 ans, pas initial de 1 an, 5 raffinements
%            directe: FickVerifMMSTemporel(10, 10, 1, 5)
%
% Historique
% 19-Oct-2022 : Création
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FickVerifMMSTemporel(Ntot, Nraf, dtInit, NdtInit)

% Validation de toutes les variables d'entrée
valid = true;
if(Ntot < 3 || (floor(Ntot) ~= ceil(Ntot)))
   disp("Le nombre de noeuds Ntot doit être un entier >= 3");
   valid = false;
end
if(dtInit <= 0)
   disp("L'intervalle de temps initial dtInit doit être > 0");
   valid = false;
end
if(NdtInit < 1 || (floor(NdtInit) ~= ceil(NdtInit)))
   disp("Le nombre de pas de temps initial NdtInit doit être une entier >= 1");
   valid = false;
end
if(Nraf < 1 || (floor(Nraf) ~= ceil(Nraf)))
   disp("Le nombre de raffinements Nraf doit être un entier >= 1");
   valid = false;
endif(~valid)
   error("Au moins un des arguments est invalide");
end

% Paramètres du problème
t0 = 0.001; % Constante temporelle [an]

% Erreurs L1, L2 et Linf
Niter = 1 + Nraf;
dt  = zeros(Niter, 1);
L1 = zeros(Niter, 1);
L2 = zeros(Niter, 1);
Linf = zeros(Niter, 2);

% Nombre de noeuds pour chaque itération
dt(1) = dtInit;
Ndt = NdtInit;
for i=1:Niter 
   if(i > 1)
      Ndt = 2*Ndt;       % Nombre de pas de temps pour cette itération
      dt(i) = dt(i-1)/2; % Pas de temps pour cette itération
   end

   % Concentrations calculées par différences finies [mol/m^3]
   [Cdf, Cmms] = FickDFMMS(Ntot, dt(i), Ndt, t0); 

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
   d = log(dt(p)/dt(p+1));
   p_hat_L1(p)   = log(L1(p)/L1(p+1))/d;
   p_hat_L2(p)   = log(L2(p)/L2(p+1))/d;
   p_hat_Linf(p) = log(Linf(p)/Linf(p+1))/d;
   disp(sprintf("pentes %d-%d: L1=%f, L2=%f, Linf=%f", p, p+1,...
        p_hat_L1(p), p_hat_L2(p), p_hat_Linf(p)));
end

% Création du graphe
figure
p1 = loglog(dt, L1, '-s');
title(sprintf("Convergence temporelle\n Ntot=%d", Ntot));
xlabel('dt [an]');
ylabel('Erreur [mol/m^{3}]');
xlim padded
ylim padded
grid on
hold on
p2 = loglog(dt, L2, '-s');
p3 = loglog(dt, Linf, '-s');
p1(1).MarkerFaceColor = p1(1).Color;
p2(1).MarkerFaceColor = p2(1).Color;
p3(1).MarkerFaceColor = p3(1).Color;

% Triangle avec les pentes des points a et b (a > b)
a = 5; % 1er point
b = 6; % 2e point
triang_x = [dt(b), dt(a)];
triang_y = exp(interp1(log(dt), log(L1), log(triang_x)));
triang_y = triang_y - [30*dt(b)^p_hat_L1(a) 30*dt(a)^p_hat_L1(a)] ;
loglog(triang_x([1,2,2,1]), triang_y([1,1,2,1]), 'k');
text(dt(a)+1e-5, exp(0.5*(log(triang_y(1))+log(triang_y(2)))), ...
     sprintf('L_{1}=%.2f\nL_{2}=%.2f\nL_{\\infty}=%.2f', ...
     p_hat_L1(a), p_hat_L2(a),p_hat_Linf(a)))

% Légende
lgd = legend('L_{l}', 'L_{2}', 'L_{\infty}');
lgd.Location = 'northwest';

end