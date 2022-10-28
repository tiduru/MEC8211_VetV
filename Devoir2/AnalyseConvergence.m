%---------------------------
% Obtention de données pour
% analyse de convergence
%---------------------------
clear all; clc; close all;
%% Settings
N = [3,5,9,17,33];
dt = [0.001,0.01,0.1,1];
dt_text = ["0001","001","01","1"];
tf = [1,10,25];

%% Data from COMSOL
DataDirectory = "./DataAnalyseConvergence/";
VS = "Variation spatiale/";
VT = "Variation temporelle/";

DataS = {}; % {N, t}(2)

for i = 1:length(N)
    for j = 1:length(tf)
        DataS{i,j} = load(DataDirectory+VS+"N"+string(N(i))+"t"+string(tf(j))+"dt"+dt_text(3)+".txt");
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
        [SimS{i,j}(:,2),SimS{i,j}(:,1)] = FickDF(N(i),dt(3),tf(j)/dt(3),2,1);
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

for i = 1:3
    figure(i);
    hold on
    loglog(0.5./(N-1), L1S(:,i));
    loglog(0.5./(N-1), L2S(:,i));
    loglog(0.5./(N-1), LinfS(:,i));
    title({"Erreur en fonction de la variation spatiale","t = "+string(tf(i))})
    xlabel("h [mm]")
    ylabel("Erreur")
    legend({"L1","L2","L_{inf}"},"location","best")
    hold off
    grid on
end

for i = 4:6
    figure(i);
    hold on
    loglog(dt, L1T(i-3,:));
    loglog(dt, L2T(i-3,:));
    loglog(dt, LinfT(i-3,:));
    title({"Erreur en fonction de la variation temporelle","h = "+string(0.5/32)+" m ; t = "+string(tf(i-3))})
    xlabel("dt [s]")
    ylabel("Erreur")
    legend({"L1","L2","L_{inf}"},"location","best")
    hold off
    grid on
end