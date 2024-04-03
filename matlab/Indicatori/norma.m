clear
close all
clc

%% Import dati

path="..\dati\curvaU_forte\";
[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

%% Impostazioni

% i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% rotazione per far coincidere il sistema di riferimento dell'accelerometro
% con quello della bicicletta
rilievo=6;


%% Import Dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db);

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

% diff_t=zeros(length(t),1);
% diff_t(1)=t(1);
% for i=2:length(t)
%     diff_t(i)=t(i)-t(i-1);
% end
% 
% disp("max diff_t: "+num2str(max(diff_t(2:end))));
% disp("min diff_t: "+num2str(min(diff_t(2:end))));


%% Accelerazione
acc=db(inizio:fine,2:4)*gzRot*9.81/-gMedio;

figure("Name","Accelerazione")
subplot(3,1,1)
plot(t,acc(:,1),LineWidth=1,Color="r")
title("Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Accelerazione Filtrata
filtered_acc=lowpass(acc,0.5,sr);

figure("Name","Accelerazione Filtrata")
subplot(3,1,1)
plot(t,filtered_acc(:,1),LineWidth=1,Color="r")
title("Accelerazione Filtrata")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Norma Accelerazione
norma_acc=zeros(length(acc),1);
norma_accXY=zeros(length(acc),1);

for i=1:length(acc)
    norma_acc(i)=norm(acc(i,:));
    norma_accXY(i)=norm(acc(i,1:2));
end

norma_filtered_acc=lowpass(norma_acc,0.5,sr);
norma_filtered_accXY=lowpass(norma_accXY,0.5,sr);


figure(Name="Norma Accelerazione")
plot(t,norma_acc,LineWidth=1,Color="r");
title("Norma Accelerazione")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Accelerazione XY")
plot(t,norma_accXY,LineWidth=1,Color="r");
title("Norma Accelerazione XY")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Accelerazione Fitrata")
plot(t,norma_filtered_acc,LineWidth=1,Color="r");
title("Norma Accelerazione Fitrata")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Accelerazione XY Filtrata")
plot(t,norma_filtered_accXY,LineWidth=1,Color="r");
title("Norma Accelerazione XY Filtrata")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Velocità
vel=cumsum(acc)*0.04;

figure("Name","Velocità")
subplot(3,1,1)
plot(t,vel(:,1),LineWidth=1,Color="r")
title("Velocità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s")
grid


%% Norma Velocità
norma_vel=zeros(length(vel),1);

for i=1:length(vel)
    norma_vel(i)=norm(vel(i,1:2));
end


figure(Name="Norma Velocità")
plot(t,norma_vel,LineWidth=1,Color="r");
title("Norma Velocità")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s")
grid


%% Velocità Angolare
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


%% Norma Velocità Angolare
norma_vang=zeros(length(vang),1);

for i=1:length(vang)
    norma_vang(i)=norm(vang(i,:));
end

norma_filtered_vang=lowpass(norma_vang,0.5,sr);


figure(Name="Norma Velocità Angolare")
plot(t,norma_vang,LineWidth=1,Color="r");
title("Norma Velocità Angolare")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Velocità Angolare Fitrata")
plot(t,norma_filtered_vang,LineWidth=1,Color="r");
title("Norma Velocità Angolare Fitrata")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Integrale Velocità Angolare
ang=cumsum(vang)*0.04;

figure("Name","Integrale Velocità Angolare")
subplot(3,1,1)
plot(t,ang(:,1),LineWidth=1,Color="r")
title("Integrale Velocità Angolare")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad")
grid


%% Norma Integrale Velocità Angolare
norma_ang=zeros(length(ang),1);
norma_angXY=zeros(length(ang),1);

for i=1:length(acc)
    norma_ang(i)=norm(ang(i,:));
    norma_angXY(i)=norm(ang(i,1:2));
end


figure(Name="Norma Integrale Velocità Angolare")
plot(t,norma_ang,LineWidth=1,Color="r");
title("Norma Integrale Velocità Angolare")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Integrale Velocità Angolare XY")
plot(t,norma_angXY,LineWidth=1,Color="r");
title("Norma Integrale Velocità Angolare XY")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Campo Magnetico
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


%% Norma Campo Magnetico
norma_mag=zeros(length(mag),1);

for i=1:length(mag)
    norma_mag(i)=norm(mag(i,:));
end

norma_filtered_mag=lowpass(norma_mag,0.5,sr);


figure(Name="Norma Campo Magnetico")
plot(t,norma_mag,LineWidth=1,Color="r");
title("Norma Campo Magnetico")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


figure(Name="Norma Campo Magnetico Fitrato")
plot(t,norma_filtered_mag,LineWidth=1,Color="r");
title("Norma Campo Magnetico Fitrato")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid

