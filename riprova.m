close all
clear
clc

%import dati accelerazioni
% dbA=importdata("dbdm\0.csv");
% acc=dbA.data(250:end-250,2:4);

dbA=importdata("dbdm\nuovo\0.csv");
acc=dbA.data(250:end-325,2:4);

%estrazione tempi
t=dbA.data(250:end-325,1)*1e-3;
t=t-t(1);


% figure
% plot3([0,acc(1,1)],[0,acc(1,2)],[0,acc(1,3)]);
% hold on
% grid
% axis equal
% plot3([0,500],[0,0],[0,0],Color="b")
% plot3([0,0],[0,500],[0,0],Color="r")
% plot3([0,0],[0,0],[0,500],Color="g")


%plot
% figure                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          igure(Name="Accelereazione")
% subplot(3,1,1);
% plot(t,acc(:,1),LineWidth=1,Color="b");
% subtitle("X", Color="b");
% subplot(3,1,2);
% plot(t,acc(:,2),LineWidth=1,Color="r");
% subtitle("Y", Color="b");
% subplot(3,1,3);
% plot(t,acc(:,3),LineWidth=1,Color="g");
% subtitle("Z", Color="b");

media=mean(acc(1:end,:));

% disp("Originale:")
% disp("  - media x: " + num2str(media(1)));
% disp("  - media y: " + num2str(media(2)));
% disp("  - media z: " + num2str(media(3)));

%% Rotazione attorno all'asse X

modYZ=sqrt(media(2)^2+media(3)^2);
thetaYZ=-acos(-media(3)/modYZ);
mRotX=[1,0,0;0,cos(thetaYZ),-sin(thetaYZ);0,sin(thetaYZ),cos(thetaYZ)];

newMedia = media*mRotX;

% disp("Ruotata attorno all'asse X:")
% disp("  - media x: " + num2str(newMedia(1)));
% disp("  - media y: " + num2str(newMedia(2)));
% disp("  - media z: " + num2str(newMedia(3)));


%% Rotazione attorno all'asse Y

modXz=sqrt(newMedia(1)^2+newMedia(3)^2);
thetaXz=-acos(-newMedia(3)/modXz);
mRoty=[cos(thetaXz),0,sin(thetaXz);0,1,0;-sin(thetaXz),0,cos(thetaXz)];

lastMedia = newMedia*mRoty;

% disp("Ruotata attorno all'asse y:")
% disp("  - media x: " + num2str(lastMedia(1)));
% disp("  - media y: " + num2str(lastMedia(2)));
% disp("  - media z: " + num2str(lastMedia(3)));

%% Matrice di rotazione
% mRot=mRotX;
mRot=mRotX*mRoty;
% mRot=mRotX*mRoty*[-1,0,0;0,-1,0;0,0,1];

%% Accelerazione di prova
g=norm(media)*1e-2;
%import dati accelerazioni
% % dbA=importdata("dbdm\BlueCoin_Log_N004.csv");
% % acc=dbA.data(2:end,2:4);

dbA=importdata("dbdm\nuovo\BlueCoin_Log_N001.csv");
acc=dbA.data(2:end,2:4)*g*1e-3;

% acc=acc.*9.81/(sqrt(sum(newMedia.^2)));

% acc=lowpass(acc,40,1e3);
% acc=highpass(acc,1,1e3);

%estrazione tempi
t=dbA.data(2:end,1)*1e-3;
t=t-t(1);

newAcc=acc*mRot;

% newAccMean=movmean(newAcc,50);
% acc=newAcc;
% newAcc=newAcc-newAccMean;

% angY=-real(acos(newAcc(:,3)/-9.81));
% for i=1:length(acc)
%     newAcc(i,:)=newAcc(i,:)*[cos(angY(i)), 0, sin(angY(i)); 0, 1, 0; -sin(angY(i)), 0, cos(angY(i))]*[-1,0,0;0,-1,0;0,0,1];
% end

newAcc=lowpass(newAcc,40,1e3);

%plot Accelerazione
figure(Name="Accelereazione")
subplot(3,2,1);
plot(t,acc(:,1),LineWidth=1,Color="b");
grid
title("Originale");
subtitle("X", Color="b");
subplot(3,2,3);
plot(t,acc(:,2),LineWidth=1,Color="r");
grid
subtitle("Y", Color="r");
subplot(3,2,5);
plot(t,acc(:,3),LineWidth=1,Color="g");
grid
subtitle("Z", Color="g");

subplot(3,2,2);
plot(t,newAcc(:,1),LineWidth=1,Color="b");
grid
title("Routato");
subtitle("X", Color="b");
subplot(3,2,4);
plot(t,newAcc(:,2),LineWidth=1,Color="r");
grid
subtitle("Y", Color="r");
subplot(3,2,6);
plot(t,newAcc(:,3),LineWidth=1,Color="g");
grid
subtitle("Z", Color="g");


%% Calcolo velocità di prova

vel=cumsum(newAcc)*0.04;

%plot
figure(Name="Velocità assi ruotati")
subplot(3,2,1);
plot(t,newAcc(:,1),LineWidth=1,Color="b");
grid
title("Accelerazione");
subtitle("X", Color="b");
subplot(3,2,3);
plot(t,newAcc(:,2),LineWidth=1,Color="r");
grid
subtitle("Y", Color="r");
subplot(3,2,5);
plot(t,newAcc(:,3),LineWidth=1,Color="g");
grid
subtitle("Z", Color="g");

subplot(3,2,2);
plot(t,vel(:,1),LineWidth=1,Color="b");
grid
title("Velocità");
subtitle("X", Color="b");
subplot(3,2,4);
plot(t,vel(:,2),LineWidth=1,Color="r");
grid
subtitle("Y", Color="r");
subplot(3,2,6);
plot(t,vel(:,3),LineWidth=1,Color="g");
grid
subtitle("Z", Color="g");


% %% Velocità angolare
% 
% gyr=dbA.data(2:end,5:7)*(1e-3*2*pi/360);
% newGyr=gyr*mRot;
% 
% % for i=1:length(acc)
% %     newGyr(i,:)=newGyr(i,:)*[cos(angY(i)), 0, sin(angY(i)); 0, 1, 0; -sin(angY(i)), 0, cos(angY(i))];
% % end
% 
% %plot Velocità Angolare
% figure(Name="Velocità Angolare")
% subplot(3,2,1);
% plot(t,gyr(:,1),LineWidth=1,Color="b");
% grid
% title("Originale");
% subtitle("X", Color="b");
% subplot(3,2,3);
% plot(t,gyr(:,2),LineWidth=1,Color="r");
% grid
% subtitle("Y", Color="r");
% subplot(3,2,5);
% plot(t,gyr(:,3),LineWidth=1,Color="g");
% grid
% subtitle("Z", Color="g");
% 
% subplot(3,2,2);
% plot(t,newGyr(:,1),LineWidth=1,Color="b");
% grid
% title("Routato");
% subtitle("X", Color="b");
% subplot(3,2,4);
% plot(t,newGyr(:,2),LineWidth=1,Color="r");
% grid
% subtitle("Y", Color="r");
% subplot(3,2,6);
% plot(t,newGyr(:,3),LineWidth=1,Color="g");
% grid
% subtitle("Z", Color="g");
% 
% %% Posizione Angolare
% 
% ang=cumtrapz(t,newGyr);
% 
% %plot
% figure(Name="Posizione Angolare")
% subplot(3,2,1);
% plot(t,newGyr(:,1),LineWidth=1,Color="b");
% grid
% title("Velocità Angolare");
% subtitle("X", Color="b");
% subplot(3,2,3);
% plot(t,newGyr(:,2),LineWidth=1,Color="r");
% grid
% subtitle("Y", Color="r");
% subplot(3,2,5);
% plot(t,newGyr(:,3),LineWidth=1,Color="g");
% grid
% subtitle("Z", Color="g");
% 
% subplot(3,2,2);
% plot(t,ang(:,1),LineWidth=1,Color="b");
% grid
% title("Posizione Angolare");
% subtitle("X", Color="b");
% subplot(3,2,4);
% plot(t,ang(:,2),LineWidth=1,Color="r");
% grid
% subtitle("Y", Color="r");
% subplot(3,2,6);
% plot(t,ang(:,3),LineWidth=1,Color="g");
% grid
% subtitle("Z", Color="g");


