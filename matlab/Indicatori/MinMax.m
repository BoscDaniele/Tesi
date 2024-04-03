clear
close all
clc

%% Import dati

path="..\dati\lunga_forte\";
[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

%% Impostazioni

% i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% rotazione per far coincidere il sistema di riferimento dell'accelerometro
% con quello della bicicletta
rilievo=4;


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


%% Massimo Accelerazione
n_acc=[10,0];
acc_max=movmax(acc,n_acc);

figure("Name","Massimo Accelerazione")
subplot(3,1,1)
plot(t,acc_max(:,1),LineWidth=1,Color="r")
title("Massimo Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Minimo Accelerazione
n_acc=[10,0];
acc_min=movmin(acc,n_acc);

figure("Name","Minimo Accelerazione")
subplot(3,1,1)
plot(t,acc_min(:,1),LineWidth=1,Color="r")
title("Minimo Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_min(:,3),LineWidth=1,Color="b")
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


%% Massimo Accelerazione Filtrata
n_fitered_acc=[10,0];
filtered_acc_max=movmax(filtered_acc,n_fitered_acc);

figure("Name","Massimo Accelerazione Filtrata")
subplot(3,1,1)
plot(t,filtered_acc_max(:,1),LineWidth=1,Color="r")
title("Massimo Accelerazione Filtrata")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Minimo Accelerazione Filtrata
n_fitered_acc=[10,0];
filtered_acc_min=movmin(filtered_acc,n_fitered_acc);

figure("Name","Minimo Accelerazione Filtrata")
subplot(3,1,1)
plot(t,filtered_acc_min(:,1),LineWidth=1,Color="r")
title("Minimo Accelerazione Filtrata")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc_min(:,3),LineWidth=1,Color="b")
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


%% Massimo Velocità
n_vel=[10,0];
vel_max=movmax(vel,n_vel);

figure("Name","Massimo Velocità")
subplot(3,1,1)
plot(t,vel_max(:,1),LineWidth=1,Color="r")
title("Massimo Velocità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s")
grid


%% Minimo Velocità
n_vel=[10,0];
vel_min=movmin(vel,n_vel);

figure("Name","Minimo Velocità")
subplot(3,1,1)
plot(t,vel_min(:,1),LineWidth=1,Color="r")
title("Minimo Velocità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel_min(:,3),LineWidth=1,Color="b")
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


%% Massimo Velocità Angolare
n_vang=[10,0];
vang_max=movmax(vang,n_vang);

figure("Name","Massimo Velocità Angolare")
subplot(3,1,1)
plot(t,vang_max(:,1),LineWidth=1,Color="r")
title("Massimo Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,vang_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,vang_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("rad/s")
grid


%% Minimo Velocità Angolare
n_vang=[10,0];
vang_min=movmin(vang,n_vang);

figure("Name","Minimo Velocità Angolare")
subplot(3,1,1)
plot(t,vang_min(:,1),LineWidth=1,Color="r")
title("Minimo Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,vang_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,vang_min(:,3),LineWidth=1,Color="b")
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


%% Massimo Integrale Velocità Angolare
n_ang=[10,0];
ang_max=movmax(ang,n_ang);

figure("Name","Massimo Integrale Velocità Angolare")
subplot(3,1,1)
plot(t,ang_max(:,1),LineWidth=1,Color="r")
title("Massimo Integrale Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ang_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ang_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("rad")
grid


%% Minimo Integrale Velocità Angolare
n_ang=[10,0];
ang_min=movmin(ang,n_ang);

figure("Name","Minimo Integrale Velocità Angolare")
subplot(3,1,1)
plot(t,ang_min(:,1),LineWidth=1,Color="r")
title("Minimo Integrale Velocità Angolare")
subtitle("X")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ang_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ang_min(:,3),LineWidth=1,Color="b")
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


%% Massimo Campo Magnetico
n_mag=[10,0];
mag_max=movmax(mag,n_mag);

figure("Name","Massimo Campo Magnetico")
subplot(3,1,1)
plot(t,mag_max(:,1),LineWidth=1,Color="r")
title("Massimo Campo Magnetico")
subtitle("X")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,2)
plot(t,mag_max(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,3)
plot(t,mag_max(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("µT")
grid


%% Minimo Campo Magnetico
n_mag=[10,0];
mag_min=movmin(mag,n_mag);

figure("Name","Minimo Campo Magnetico")
subplot(3,1,1)
plot(t,mag_min(:,1),LineWidth=1,Color="r")
title("Minimo Campo Magnetico")
subtitle("X")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,2)
plot(t,mag_min(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,3)
plot(t,mag_min(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("µT")
grid


