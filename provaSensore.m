close all
clear
clc

%import dati accelerazioni
dbA=importdata("dbdm\BlueCoin_Log_N003.csv");
acc=dbA.data(2:end,2:4);

%estrazione tempi
t=dbA.data(2:end,1);
t=t-t(1);

%import dati rotazioni
dbR=importdata("dbdm\BlueCoin_Log_N003.csv");
rot=dbR.data(2:end,5:7);


%calcolo matrice di rotazione
[mRot,mRotY,mRotx]=RotMat(acc(1,:));


%% Accelerazioni

% calcolo nuova matrice delle accelerazioni
newAcc = zeros(length(acc),3);
for i=1:length(acc)
    newAcc(i,:)=acc(i,:)*mRot;
end

%grafico 3D vettore gravità
figure
plot3([0,acc(1,1)],[0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="black")
grid
hold on
% plot3([0,newAcc(1,1)],[0,newAcc(1,2)],[0,newAcc(1,3)],LineWidth=1,Color="black")
plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="b")
plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="r")
plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="g")
xlabel("X", Color="b")
ylabel("Y",Color="r")
zlabel("Z", Color="g")

disp("modulo = " + sprintf(num2str(sqrt(acc(1,1)^2 + acc(1,2)^2 + acc(1,3)^2))));

%grafici accelerazioni
figure(Name="Accelerazioni")
subplot(3,2,1);
plot(t,acc(:,1),LineWidth=1,Color="b");
title("Originale");
subtitle("X", Color="b");

subplot(3,2,2);
plot(t,newAcc(:,1),LineWidth=1,Color="b");
title("MenoOriginale");
subtitle("X", Color="b");

subplot(3,2,3);
plot(t,acc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,4);
plot(t,newAcc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,5);
plot(t,acc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

subplot(3,2,6);
plot(t,newAcc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

moduloPrima=zeros(length(acc),1);
moduloDopo=zeros(length(acc),1);

for i=1:length(acc)
    moduloPrima(i)=sqrt(acc(i,1)^2+acc(i,2)^2+acc(i,3)^2);
    moduloDopo(i)=sqrt(newAcc(i,1)^2+newAcc(i,2)^2+newAcc(i,3)^2);
end

figure(Name="Modulo")
subplot(1,2,1)
plot(t,moduloPrima,LineWidth=1)
subtitle("modulo prima")

subplot(1,2,2)
plot(t,moduloDopo,LineWidth=1)
subtitle("modulo dopo")

figure(Name="Differenza Moduli")
plot(t,moduloPrima-moduloDopo,LineWidth=1)
subtitle("Modulo Diff")

%% Rotazioni

%calcolo nuova matrice delle rotazioni
newRot = zeros(length(rot),3);
for i=1:length(rot)
    newRot(i,:)=rot(i,:)*mRot;
end

%grafici rotazioni
figure(Name="Rotazioni")
subplot(3,2,1);
plot(t,rot(:,1),LineWidth=1,Color="b");
title("Originale");
subtitle("X", Color="b");

subplot(3,2,2);
plot(t,-newRot(:,1),LineWidth=1,Color="b");
title("MenoOriginale");
subtitle("X", Color="b");

subplot(3,2,3);
plot(t,rot(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,4);
plot(t,-newRot(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,5);
plot(t,rot(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

subplot(3,2,6);
plot(t,-newRot(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

moduloPrimaR=zeros(length(rot),1);
moduloDopoR=zeros(length(rot),1);

for i=1:length(acc)
    moduloPrimaR(i)=sqrt(rot(i,1)^2+rot(i,2)^2+rot(i,3)^2);
    moduloDopoR(i)=sqrt(newRot(i,1)^2+newRot(i,2)^2+newRot(i,3)^2);
end

figure(Name="Modulo Rotazione")
subplot(1,2,1)
plot(t,moduloPrimaR,LineWidth=1)
subtitle("modulo prima")

subplot(1,2,2)
plot(t,moduloDopoR,LineWidth=1)
subtitle("modulo dopo")

figure(Name="Differenza Moduli Rotazione")
plot(t,moduloPrimaR-moduloDopoR,LineWidth=1)
subtitle("Modulo Diff")

