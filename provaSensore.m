close all
clear
clc

dbA=importdata("dbdm\Acceleration.csv");
acc=dbA.data;

t=str2num(cell2mat(dbA.textdata(4:end,4)));
t=t-t(1);

%Grafico componente XZ e YZ del primo vettore
% figure(Name="Primo Vettore XZ-YZ")
% subplot(1,2,1)
% plot([0,acc(1,1)],[0,acc(1,3)],LineWidth=1,Color="black")
% grid
% xlabel("x",Color="b")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,acc(1,3)],LineWidth=1,Color="g")
% plot([0,acc(1,1)],[0,0],LineWidth=1,Color="b")
% subtitle("Componente XZ",Color="b")
%
% subplot(1,2,2)
% plot([0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="black")
% grid
% xlabel("y",Color="r")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,acc(1,3)],LineWidth=1,Color="g")
% plot([0,acc(1,2)],[0,0],LineWidth=1,Color="r")
% subtitle("Componente YZ",Color="r")

%rotazione attorno all'asse Y
xz=sqrt(acc(1,1)^2+acc(1,3)^2);
thetaXZ=acos(acc(1,3)/xz);
mRotY=[cos(thetaXZ), 0, sin(thetaXZ); 0, 1, 0; -sin(thetaXZ), 0, cos(thetaXZ)];
accY=acc(1,:)*mRotY;

%Grafico componente xz e Yz del vettore ruotato attorno all'asse Y
% figure(Name="Assi routati attorno a Y")
% subplot(1,2,1)
% plot([0,accY(1)],[0,accY(3)],LineWidth=1,Color="black")
% grid
% xlabel("x",Color="b")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,accY(3)],LineWidth=1,Color="g")
% plot([0,accY(1)],[0,0],LineWidth=1,Color="b")
% subtitle("Componente xz",Color="b")
%
% subplot(1,2,2)
% plot([0,accY(2)],[0,accY(3)],LineWidth=1,Color="black")
% grid
% xlabel("y",Color="r")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,accY(3)],LineWidth=1,Color="g")
% plot([0,accY(2)],[0,0],LineWidth=1,Color="r")
% subtitle("Componente Yz",Color="r")

%rotazione attorno all'asse x
yz=sqrt(accY(2)^2+accY(3)^2);
thetaYZ=acos(accY(3)/yz);
mRotx=[1, 0, 0; 0, cos(thetaYZ), sin(thetaYZ); 0, -sin(thetaYZ), cos(thetaYZ)];
accx=accY*mRotx;

mRot=mRotY*mRotx;

%Grafico componente xz e yz del vettore ruotato attorno all'asse x
% figure(Name="Assi routati attorno a x")
% subplot(1,2,1)
% plot([0,accx(1)],[0,accx(3)],LineWidth=1,Color="black")
% grid
% xlabel("x",Color="b")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,accx(3)],LineWidth=1,Color="g")
% plot([0,accx(1)],[0,0],LineWidth=1,Color="b")
% subtitle("Componente xz",Color="b")
%
% subplot(1,2,2)
% plot([0,accx(2)],[0,accx(3)],LineWidth=1,Color="black")
% grid
% xlabel("y",Color="r")
% ylabel("z",Color="g")
% hold on
% plot([0,0],[0,accx(3)],LineWidth=1,Color="g")
% plot([0,accx(2)],[0,0],LineWidth=1,Color="r")
% subtitle("Componente yz",Color="r")

%Grafico 3D primo vettore
% figure(Name="Primo Vettore 3D")
% plot3([0,acc(1,1)],[0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="black");
% xlabel("x",Color="b")
% ylabel("y",Color="r")
% zlabel("z",Color="g")
% hold on
% grid
% plot3([0,acc(1,1)],[0,0],[0,0],LineWidth=1,Color="b");
% plot3([0,0],[0,acc(1,2)],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[0,0],[0,acc(1,3)],LineWidth=1,Color="g");

%Grafico 3D assi routati
% figure(Name="Assi Ruotati 3D")
% plot3([0,accx(1)],[0,accx(2)],[0,accx(3)],LineWidth=1,Color="black");
% xlabel("x",Color="b")
% ylabel("y",Color="r")
% zlabel("z",Color="g")
% hold on
% grid
% plot3([0,1000],[0,0],[0,0],LineWidth=1,Color="b");
% plot3([0,0],[0,1000],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[0,0],[0,1000],LineWidth=1,Color="g");

%grafico direzione primo vettore
% figure
% subplot(3,2,1)
% plot([0,acc(1,1)],[0,acc(1,2)],LineWidth=1,Color="b")
% grid
% title("Originale")
% subtitle("XY",Color="b")
%
% subplot(3,2,2)
% plot([0,accx(1)],[0,accx(2)],LineWidth=1,Color="b")
% grid
% title("Routato")
% subtitle("newXY",Color="b")
%
% subplot(3,2,3)
% plot([0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="r")
% grid
% subtitle("YZ",Color="r")
%
% subplot(3,2,4)
% plot([0,accx(2)],[0,accx(3)],LineWidth=1,Color="r")
% grid
% subtitle("newYZ",Color="r")
%
% subplot(3,2,5)
% plot([0,acc(1,1)],[0,acc(1,3)],LineWidth=1,Color="g")
% grid
% subtitle("XZ",Color="g")
%
% subplot(3,2,6)
% plot([0,accx(1)],[0,accx(3)],LineWidth=1,Color="g")
% grid
% subtitle("newXZ",Color="g")

%calcolo nuova matrice delle accelerazioni
newAcc = zeros(length(acc),3);

for i=1:length(acc)
    newAcc(i,:)=acc(i,:)*mRot;
end

%grafici accelerazioni
subplot(3,2,1);
plot(t,acc(:,1),LineWidth=1,Color="b");
title("Originale");
subtitle("X", Color="b");

subplot(3,2,2);
plot(t,-newAcc(:,1),LineWidth=1,Color="b");
title("MenoOriginale");
subtitle("X", Color="b");

subplot(3,2,3);
plot(t,acc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,4);
plot(t,-newAcc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,5);
plot(t,acc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

subplot(3,2,6);
plot(t,-newAcc(:,3),LineWidth=1,Color="g");
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