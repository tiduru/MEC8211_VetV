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

DataS{1,1} = load(DataDirectory+VS+"N"+string(N(1))+"t"+string(tf(1))+"dt"+dt_text(4)+".txt");
DataS{1,2} = load(DataDirectory+VS+"N"+string(N(1))+"t"+string(tf(2))+"dt"+dt_text(4)+".txt");
DataS{1,3} = load(DataDirectory+VS+"N"+string(N(1))+"t"+string(tf(3))+"dt"+dt_text(4)+".txt");

DataS{2,1} = load(DataDirectory+VS+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(4)+".txt");
DataS{2,2} = load(DataDirectory+VS+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(4)+".txt");
DataS{2,3} = load(DataDirectory+VS+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(4)+".txt");

DataS{3,1} = load(DataDirectory+VS+"N"+string(N(3))+"t"+string(tf(1))+"dt"+dt_text(4)+".txt");
DataS{3,2} = load(DataDirectory+VS+"N"+string(N(3))+"t"+string(tf(2))+"dt"+dt_text(4)+".txt");
DataS{3,3} = load(DataDirectory+VS+"N"+string(N(3))+"t"+string(tf(3))+"dt"+dt_text(4)+".txt");

DataS{4,1} = load(DataDirectory+VS+"N"+string(N(4))+"t"+string(tf(1))+"dt"+dt_text(4)+".txt");
DataS{4,2} = load(DataDirectory+VS+"N"+string(N(4))+"t"+string(tf(2))+"dt"+dt_text(4)+".txt");
DataS{4,3} = load(DataDirectory+VS+"N"+string(N(4))+"t"+string(tf(3))+"dt"+dt_text(4)+".txt");

DataS{5,1} = load(DataDirectory+VS+"N"+string(N(5))+"t"+string(tf(1))+"dt"+dt_text(4)+".txt");
DataS{5,2} = load(DataDirectory+VS+"N"+string(N(5))+"t"+string(tf(2))+"dt"+dt_text(4)+".txt");
DataS{5,3} = load(DataDirectory+VS+"N"+string(N(5))+"t"+string(tf(3))+"dt"+dt_text(4)+".txt");

DataT5 = {}; % {t, dt}(2)

DataT5{1,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(1)+".txt");
DataT5{1,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(2)+".txt");
DataT5{1,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(3)+".txt");
DataT5{1,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(5)+".txt");

DataT5{2,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(1)+".txt");
DataT5{2,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(2)+".txt");
DataT5{2,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(3)+".txt");
DataT5{2,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(5)+".txt");

DataT5{3,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(1)+".txt");
DataT5{3,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(2)+".txt");
DataT5{3,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(3)+".txt");
DataT5{3,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(5)+".txt");

DataT17 = {}; % {t,dt}(2)

DataT17{1,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(1)+".txt");
DataT17{1,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(2)+".txt");
DataT17{1,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(3)+".txt");
DataT17{1,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(1))+"dt"+dt_text(5)+".txt");

DataT17{2,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(1)+".txt");
DataT17{2,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(2)+".txt");
DataT17{2,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(3)+".txt");
DataT17{2,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(2))+"dt"+dt_text(5)+".txt");

DataT17{3,1} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(1)+".txt");
DataT17{3,2} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(2)+".txt");
DataT17{3,3} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(3)+".txt");
DataT17{3,4} = load(DataDirectory+VT(1)+"N"+string(N(2))+"t"+string(tf(3))+"dt"+dt_text(5)+".txt");

%% Settings for own simulations

