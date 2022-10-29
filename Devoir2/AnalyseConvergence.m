%---------------------------
% Obtention de données pour
% analyse de convergence
%---------------------------
clear all; clc; close all;
%% Settings
N = [3,5,9,17,33];
h = 0.5./(N-1);
dt = [0.001,0.01,0.1,1];
dt_text = ["0001","001","01","1"];
tf = [10,50];

%% Data from COMSOL
DataDirectory = "./DataAnalyseConvergence/";
VS = "Variation spatiale/";
VT = "Variation temporelle/";

DataS = {}; % {N, t}(2)

for i = 1:length(N)
    for j = 1:length(tf)
        DataS{i,j} = load(DataDirectory+VS+"N"+string(N(i))+"t"+string(tf(j))+"dt"+dt_text(1)+".txt");
    end
end

DataT = {}; % {t, dt}(2)

for i = 1:length(tf)
    for j = 1:length(dt)
        DataT{i,j} = load(DataDirectory+VT+"N"+string(N(end))+"t"+string(tf(i))+"dt"+dt_text(j)+".txt");
    end
end

%% Settings for own simulations
SimS = {};
for i = 1:length(N)
    for j = 1:length(tf)
        [SimS{i,j}(:,2),SimS{i,j}(:,1)] = FickDF(N(i),dt(1),tf(j)/dt(1),2,1);
    end
end

SimT = {};
for i = 1:length(tf)
    for j = 1:length(dt)
        [SimT{i,j}(:,2) , SimT{i,j}(:,1)] = FickDF(N(end),dt(j),tf(i)/dt(j),2,1);
    end
end

%% Error calculations
L1S = [];
L2S = [];
LinfS = [];

L1T = [];
L2T = [];
LinfT = [];

for i = 1:length(N)
    for j = 1:length(tf)
        L1S(i,j) = 0;
        L2S(i,j) = 0;
        LinfS(i,j) = 0;
        for k = 1:size(DataS{i,j},2)
            diffabs = abs(DataS{i,j}(k,2) - SimS{i,j}(k,2));
            L1S(i,j) = L1S(i,j) + diffabs;
            L2S(i,j) = L2S(i,j) + diffabs^2;
            if diffabs > LinfS(i,j)
                LinfS(i,j) = diffabs;
            end
        end
        L1S(i,j) = L1S(i,j) / N(i);
        L2S(i,j) = sqrt(L2S(i,j) / N(i));
    end
end

for i = 1:length(tf)
    for j = 1:length(dt)
        L1T(i,j) = 0;
        L2T(i,j) = 0;
        LinfT(i,j) = 0;
        for k = 1:size(DataT{i,j},2)
            diffabs = abs(DataT{i,j}(k,2) - SimT{i,j}(k,2));
            L1T(i,j) = L1T(i,j) + diffabs;
            L2T(i,j) = L2T(i,j) + diffabs^2;
            if diffabs > LinfT(i,j)
                LinfT(i,j) = diffabs;
            end
        end
        L1T(i,j) = L1T(i,j) / N(2);
        L2T(i,j) = sqrt(L2T(i,j) / N(2));
    end
end

% Slopes
p_hat_L1S = [];
p_hat_L2S = [];
p_hat_LinfS = [];

p_hat_L1T = [];
p_hat_L2T = [];
p_hat_LinfT = [];

disp("**Convergence spatiale**");
for j = 1:length(tf)
   for i = 1:length(N)-1
       dh = log(h(i)/h(i+1));
       p_hat_L1S(i,j)   = log(L1S(i,j)/L1S(i+1,j))/dh;
       p_hat_L2S(i,j)   = log(L2S(i,j)/L2S(i+1,j))/dh;
       p_hat_LinfS(i,j) = log(LinfS(i,j)/LinfS(i+1,j))/dh;
       disp(sprintf("Time=%.3f, Slopes %d-%d: L1=%.2f, L2=%.2f, Linf=%.2f", ...
            tf(j), i, i+1, p_hat_L1S(i,j), p_hat_L2S(i,j), p_hat_LinfS(i,j)));
    end
end

disp("**Convergence temporelle**");
for i = 1:length(tf)
    for j = 1:length(dt)-1
       ddt = log(dt(j)/dt(j+1));
       p_hat_L1T(i,j)   = log(L1T(i,j)/L1T(i,j+1))/ddt;
       p_hat_L2T(i,j)   = log(L2T(i,j)/L2T(i,j+1))/ddt;
       p_hat_LinfT(i,j) = log(LinfT(i,j)/LinfT(i,j+1))/ddt;
       disp(sprintf("Time=%.3f, Slopes %d-%d: L1=%.2f, L2=%.2f, Linf=%.2f", ...
            tf(i), j, j+1, p_hat_L1T(i,j), p_hat_L2T(i,j), p_hat_LinfT(i,j)));
    end
end

for i = 1:2
    figure(i);
    axes("XScale","log","YScale","log");
    set(gca,'Box','on');
    hold on
    p1 = loglog(h, L1S(:,i), '-s');
    p2 = loglog(h, L2S(:,i), '-s');
    p3 = loglog(h, LinfS(:,i), '-s');
    
    p1(1).MarkerFaceColor = p1(1).Color;
    p2(1).MarkerFaceColor = p2(1).Color;
    p3(1).MarkerFaceColor = p3(1).Color;
    
    title({"Erreur en fonction de la variation spatiale","t = "+string(tf(i))})
    xlabel("h [m]")
    ylabel("Erreur [mol/m^{3}]")
    xlim padded
    ylim padded
    grid on

    % Triangle avec les pentes des points a et b (a > b)
    a = 2; % 1er point
    b = 3; % 2e point
    triang_x = [h(b), h(a)];
    triang_y = exp(interp1(log(h), log(L1S(:,i)), log(triang_x)));
    triang_y = triang_y - 5*[h(b)^p_hat_L1S(a,i) h(a)^p_hat_L1S(a,i)] ;
    loglog(triang_x([1,2,2,1]), triang_y([1,1,2,1]), 'k');
    text(h(a)+0.001, exp(0.5*(log(triang_y(1))+log(triang_y(2)))), ...
    sprintf('L_{1}=%.2f\nL_{2}=%.2f\nL_{\\infty}=%.2f', ...
    p_hat_L1S(a,i), p_hat_L2S(a,i),p_hat_LinfS(a,i)))

    lgd = legend("L1","L2","L_{\infty}", "");
    lgd.Location = 'northwest';

    hold off
end

for i = 3:4
    figure(i);
    axes("XScale","log","YScale","log");
    set(gca,'Box','on');
    hold on
    p1 = loglog(dt, L1T(i-2,:));
    p2 = loglog(dt, L2T(i-2,:));
    p3 = loglog(dt, LinfT(i-2,:));

    p1(1).MarkerFaceColor = p1(1).Color;
    p2(1).MarkerFaceColor = p2(1).Color;
    p3(1).MarkerFaceColor = p3(1).Color;

    title({"Erreur en fonction de la variation temporelle","h = "+string(0.5/32)+" m ; t = "+string(tf(i-2))})
    xlabel("dt [an]")
    ylabel("Erreur [mol/m^{3}]")
    xlim padded
    ylim padded
    grid on
    legend({"L_{1}","L_{2}","L_{\infty}"},"location","best")
    hold off
    grid on
end