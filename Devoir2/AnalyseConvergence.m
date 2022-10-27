%---------------------------
% Obtention de données pour
% analyse de convergence
%---------------------------
clear all; clc;
%% Settings
N = [3,5,9,17,33];
dt = [0.001,0.01,0.1,0.25,1];
dt_text = ["0001","001","01","025","1"];
tf = [1,10,25];

%% Data from COMSOL
DataDirectory = "./DataAnalyseConvergence/";
VS = "Variation spatial/";
VT = ["Variation temporel/N5/","Variation temporel/N17/"];

DataS = {}; % {N, t}(2)

for i = 1:length(N)
    for j = 1:length(tf)
        DataS{i,j} = load(DataDirectory+VS+"N"+string(N(i))+"t"+string(tf(j))+"dt"+dt_text(4)+".txt");
    end
end

DataT5 = {}; % {t, dt}(2)
DataT17 = {};

for i = 1:length(tf)
    for j = 1:length(dt)
        DataT5{i,j} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(i))+"dt"+dt_text(j)+".txt");
        DataT17{i,j} = load(DataDirectory+VT(2)+"N"+string(N(4))+"t"+string(tf(i))+"dt"+dt_text(j)+".txt");
    end
end

%% Settings for own simulations
SimS = {};
for i = 1:length(N)
    for j = 1:length(tf)
        [SimS{i,j}(:,2),SimS{i,j}(:,1)] = FickDF(N(i),dt(4),tf(j)/dt(4),2,1);
    end
end

SimT5 = {};
SimT17 = {};
for i = 1:length(tf)
    for j = 1:length(dt)
        [SimT5{i,j}(:,2) , SimT5{i,j}(:,1)] = FickDF(N(2),dt(j),tf(i)/dt(j),2,1);
        [SimT17{i,j}(:,2) , SimT17{i,j}(:,1)] = FickDF(N(4),dt(j),tf(i)/dt(j),2,1);
    end
end

%% Error calculations
L1S = [];
L2S = [];
LinfS = [];

L1T5 = [];
L2T5 = [];
LinfT5 = [];

L1T17 = [];
L2T17 = [];
LinfT17 = [];

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
        L1T5(i,j) = 0;
        L2T5(i,j) = 0;
        LinfT5(i,j) = 0;
        for k = 1:size(DataT5{i,j},2)
            diffabs = abs(DataT5{i,j}(k,2) - SimT5{i,j}(k,2));
            L1T5(i,j) = L1T5(i,j) + diffabs;
            L2T5(i,j) = L2T5(i,j) + diffabs^2;
            if diffabs > LinfT5(i,j)
                LinfT5(i,j) = diffabs;
            end
        end
        L1T5(i,j) = L1T5(i,j) / N(2);
        L2T5(i,j) = sqrt(L2T5(i,j) / N(2));
    end
end
for i = 1:length(tf)
    for j = 1:length(dt)
        L1T17(i,j) = 0;
        L2T17(i,j) = 0;
        LinfT17(i,j) = 0;
        for k = 1:size(DataT17{i,j},2)
            diffabs = abs(DataT17{i,j}(k,2) - SimT17{i,j}(k,2));
            L1T17(i,j) = L1T17(i,j) + diffabs;
            L2T17(i,j) = L2T17(i,j) + diffabs^2;
            if diffabs > LinfT17(i,j)
                LinfT17(i,j) = diffabs;
            end
        end
        L1T17(i,j) = L1T17(i,j) / N(4);
        L2T17(i,j) = sqrt(L2T17(i,j) / N(4));
    end
end

for i = 1:3
    figure(i);
    hold on
    loglog(0.5./N, L1S(:,i));
    loglog(0.5./N, L2S(:,i));
    loglog(0.5./N, LinfS(:,i));
    hold off
end

for i = 4:6
    figure(i);
    hold on
    loglog(dt, L1T5(i-3,:));
    loglog(dt, L2T5(i-3,:));
    loglog(dt, LinfT5(i-3,:));
    hold off
end
