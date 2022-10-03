%% Clear l'espace de travail
clc; clear all; close all;

%%
%Déclaration des constantes
k=4*10^-9;Deff=10^(-10);S=10^-8;Ce=10;

%Création du maillage en espace et en temps
%En espace
r_A=0;r_B=0.5;
NELM_espace=5;
dr=(r_B-r_A)/(NELM_espace-1);

%En temps
T=1000;
NELM_temps=10000;
dt=T/NELM_temps;

%Création des variables à indexer
syms i; 
syms r(i); r(i)=(i-1)*dr;
syms A(i); A(i)=r*dt*Deff;
syms B(i); B(i)=r*dr^2;
syms F(i);F(i)=-S*dt*r*dr^2;
E=Deff*dt*dr;

%Initialisation des matrices permettant la résolution du système
%matricielle linéaire
A2=zeros(NELM_espace,NELM_espace);
A2(NELM_espace,NELM_espace)=1;
Solution=zeros(NELM_espace,T/dt); Solution(NELM_espace,:)=Ce;
A2(1,1:2)=[1 -1];
B2 = zeros(NELM_espace,1); B2(end,1)=Ce;

%% Résolution S constant

syms K(i) [1 3]; K=[-A(i) B(i)+E+2*A(i) -E-A(i)];
for j=2:NELM_espace-1
    B2(j,1)=F(j)+r(j)*dr^2*Solution(j,1);
    A2(j,j-1:j+1)=subs(K,i,j);
end

%Résolution dans le temps du problème par différence finis
for j=2:T/dt-1
Solution(:,j)=A2\B2;
B2(2:NELM_espace-1,1)=F(2:NELM_espace-1)+r(2:NELM_espace-1)*dr^2*Solution(2:NELM_espace-1,j);
end
Solution(:,T/dt)=A2\B2;

%% Résolution S = kC

syms K(i) [1 3]; K=[-A(i) B(i)+E+2*A(i)+k -E-A(i)];
for j=2:NELM_espace-1
    B2(j,1)=r(j)*dr^2*Solution(j,1);
    A2(j,j-1:j+1)=subs(K,i,j);
end

%Résolution dans le temps du problème par différence finis
for j=2:T/dt-1
Solution(:,j)=A2\B2;
B2(2:NELM_espace-1,1)=r(2:NELM_espace-1)*dr^2*Solution(2:NELM_espace-1,j);
end
Solution(:,T/dt)=A2\B2;

%% Plot de la solution stationnaire
syms Sol_stationnaire(i);
Sol_stationnaire(i)=S/(4*Deff)*r_B^2*(r(i)^2/r_B^2-1)+Ce;

hold on;
figure (1);
plot(subs(r,i,1:NELM_espace),subs(Sol_stationnaire,i,1:NELM_espace));
plot(subs(r,i,1:NELM_espace),Solution(:,end));