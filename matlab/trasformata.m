clear
close all
clc

%% Import dati

path="db\secchia\";
[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

%% Impostazioni

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - discesa via secchia
% 3 - salita mi fermo in salita
% 4 - discesa parto in discesa
% 5 - curve
% 6 - pedalata tranquilla
% 7 - pedalata forte
rilievo=5;


%% Import Dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db);

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*9.81/-gMedio;

acc_rot=acc*gzRot;

figure("Name","Accelerazione")
subplot(3,1,1)
plot(t,acc_rot(:,1),LineWidth=1,Color="r")
title("Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_rot(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_rot(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


vang=deg2rad(db(inizio:fine,5:7)*1e-3);

figure("Name","Velocità Angolare")
subplot(3,1,1)
plot(t,vang(:,1),LineWidth=1,Color="r")
title("Velocità Angolare")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,vang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,vang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad/s")
grid


% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);

figure("Name","Campo Magnetico")
subplot(3,1,1)
plot(t,mag(:,1),LineWidth=1,Color="r")
title("Campo Magnetico")
subtitle("X")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,2)
plot(t,mag(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,3)
plot(t,mag(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("µT")
grid


%% Trasformata Accelerazione
inizio_f=2;

L=length(acc(:,1));
f=sr/L*(0:(L/2));

Y=fft(acc);
P2=abs(Y/L);
trasform_acc=P2(1:(L/2+1),:);
trasform_acc(2:end-1,:)=2*trasform_acc(2:end-1,:);

figure("Name","Trasformata Accelerazione")
subplot(3,1,1)
plot(f(inizio_f:end),trasform_acc(inizio_f:end,1),LineWidth=1,Color="r");
title("Trasformata Accelerazione")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,1])
grid
subplot(3,1,2)
plot(f(inizio_f:end),trasform_acc(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,1])
grid
subplot(3,1,3)
plot(f(inizio_f:end),trasform_acc(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,1])
grid


%% Spettro Accelerazione
xdft=fft(acc);
xdft=xdft(1:L/2+1,:);
dens_acc=(1/(sr*L))*abs(xdft).^2;
dens_acc(2:end-1,:)=2*dens_acc(2:end-1,:);

figure("Name","Spettro Accelerazione")
subplot(3,1,1)
plot(f(inizio_f:end),dens_acc(inizio_f:end,1),LineWidth=1,Color="r");
title("Spettro Accelerazione")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,8])
grid
subplot(3,1,2)
plot(f(inizio_f:end),dens_acc(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,2])
grid
subplot(3,1,3)
plot(f(inizio_f:end),dens_acc(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,8])
grid


%% Trasformata Velocità Angolare
Y=fft(vang);
P2=abs(Y/L);
trasform_vang=P2(1:(L/2+1),:);
trasform_vang(2:end-1,:)=2*trasform_vang(2:end-1,:);

figure("Name","Trasformata Velocità Angolare")
subplot(3,1,1)
plot(f(inizio_f:end),trasform_vang(inizio_f:end,1),LineWidth=1,Color="r");
title("Trasformata Velocità Angolare")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,1])
grid
subplot(3,1,2)
plot(f(inizio_f:end),trasform_vang(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,1])
grid
subplot(3,1,3)
plot(f(inizio_f:end),trasform_vang(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,1])
grid


%% Spettro Velocità Angolare
xdft=fft(vang);
xdft=xdft(1:L/2+1,:);
dens_vang=(1/(sr*L))*abs(xdft).^2;
dens_vang(2:end-1,:)=2*dens_vang(2:end-1,:);

figure("Name","Spettro Velocità Angolare")
subplot(3,1,1)
plot(f(inizio_f:end),dens_vang(inizio_f:end,1),LineWidth=1,Color="r");
title("Spettro Spettro Velocità Angolare")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,8])
grid
subplot(3,1,2)
plot(f(inizio_f:end),dens_vang(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,2])
grid
subplot(3,1,3)
plot(f(inizio_f:end),dens_vang(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,8])
grid

