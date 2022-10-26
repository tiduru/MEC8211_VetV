%---------------------------
% Obtention de données pour
% analyse de convergence
%---------------------------

%% Data capture from COMSOL


%% Settings for own simulations
N = 3:2:33;
dt = 0.25;
tf = 25;

C = {};
r = {};

for i = 1:length(N)
    [C{i}, r{i}] = FickDF(N(i), dt, tf/dt, 2, 1);
end

