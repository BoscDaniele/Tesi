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


%% Distanza Picco Picco Accelerazione
n_acc=[10,0];
acc_peak=movmax(acc,n_acc)-movmin(acc,n_acc);

figure("Name","Distanza Picco Picco Accelerazione")
subplot(3,1,1)
plot(t,acc_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_peak(:,3),LineWidth=1,Color="b")
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


%% Distanza Picco Picco Accelerazione Filtrata
n_fitered_acc=[10,0];
filtered_acc_peak=movmax(filtered_acc,n_fitered_acc)-movmin(filtered_acc,n_fitered_acc);

figure("Name","Distanza Picco Picco Accelerazione Filtrata")
subplot(3,1,1)
plot(t,filtered_acc_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Accelerazione Filtrata")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc_peak(:,3),LineWidth=1,Color="b")
subtitle("Z")
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


%% Distanza Picco Picco Velocità
n_vel=[10,0];
vel_peak=movmax(vel,n_vel)-movmin(vel,n_vel);

figure("Name","Distanza Picco Picco Velocità")
subplot(3,1,1)
plot(t,vel_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Velocità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel_peak(:,3),LineWidth=1,Color="b")
subtitle("Z")
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


%% Distanza Picco Picco Velocità Angolare
n_vang=[10,0];
vang_peak=movmax(vang,n_vang)-movmin(vang,n_vang);

figure("Name","Distanza Picco Picco Velocità Angolare")
subplot(3,1,1)
plot(t,vang_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,vang_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,vang_peak(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("rad/s")
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


%% Distanza Picco Picco Integrale Velocità Angolare
n_ang=[10,0];
ang_peak=movmax(ang,n_ang)-movmin(ang,n_ang);

figure("Name","Distanza Picco Picco Integrale Velocità Angolare")
subplot(3,1,1)
plot(t,ang_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Integrale Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ang_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ang_peak(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("rad")
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


%% Distanza Picco Picco Campo Magnetico
n_mag=[10,0];
mag_peak=movmax(mag,n_mag)-movmin(mag,n_mag);

figure("Name","Distanza Picco Picco Campo Magnetico")
subplot(3,1,1)
plot(t,mag_peak(:,1),LineWidth=1,Color="r")
title("Distanza Picco Picco Campo Magnetico")
subtitle("X")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,2)
plot(t,mag_peak(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,3)
plot(t,mag_peak(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("µT")
grid


