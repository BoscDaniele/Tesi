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


%% Densità Spettrale Accelerazione
inizio_f=1;

L=length(acc(:,1));
freq=0:2*pi/L:pi;

xdft=fft(acc);
xdft=xdft(1:L/2+1,:);
dens_acc=(1/(2*pi*L))*abs(xdft).^2;
dens_acc(2:end-1,:)=2*dens_acc(2:end-1,:);

figure("Name","Densità Spettrale di Potenza Accelerazione")
subplot(3,1,1)
plot(freq,pow2db(dens_acc(:,1)),LineWidth=1,Color="r");
title("Densità Spettrale di Potenza Accelerazione")
subtitle("X''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid
subplot(3,1,2)
plot(freq,pow2db(dens_acc(:,2)),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,2])
grid
subplot(3,1,3)
plot(freq,pow2db(dens_acc(:,3)),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
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


%% Densità Spettrale Velocità
xdft=fft(vel);
xdft=xdft(1:L/2+1,:);
dens_vel=(1/(2*pi*L))*abs(xdft).^2;
dens_vel(2:end-1,:)=2*dens_vel(2:end-1,:);

figure("Name","Densità Spettrale Velocità")
subplot(3,1,1)
plot(freq,pow2db(dens_vel(:,1)),LineWidth=1,Color="r");
title("Densità Spettrale Velocità")
subtitle("X''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid
subplot(3,1,2)
plot(freq,pow2db(dens_vel(:,2)),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,2])
grid
subplot(3,1,3)
plot(freq,pow2db(dens_vel(:,3)),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
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


%% Densità Spettrale Velocità Angolare
xdft=fft(vang);
xdft=xdft(1:L/2+1,:);
dens_vang=(1/(2*pi*L))*abs(xdft).^2;
dens_vang(2:end-1,:)=2*dens_vang(2:end-1,:);

figure("Name","Densità Spettrale Velocità Angolare")
subplot(3,1,1)
plot(freq,pow2db(dens_vang(:,1)),LineWidth=1,Color="r");
title("Densità Spettrale Velocità Angolare")
subtitle("X''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid
subplot(3,1,2)
plot(freq,pow2db(dens_vang(:,2)),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,2])
grid
subplot(3,1,3)
plot(freq,pow2db(dens_vang(:,3)),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
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


%% Spettro Integrale Velocità Angolare
xdft=fft(ang);
xdft=xdft(1:L/2+1,:);
dens_ang=(1/(2*pi*L))*abs(xdft).^2;
dens_ang(2:end-1,:)=2*dens_ang(2:end-1,:);

figure("Name","Densità Spettrale Integrale Veloctià Angolare")
subplot(3,1,1)
plot(freq,pow2db(dens_ang(:,1)),LineWidth=1,Color="r");
title("Densità Spettrale Velocità Angolare")
subtitle("X''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid
subplot(3,1,2)
plot(freq,pow2db(dens_ang(:,2)),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,2])
grid
subplot(3,1,3)
plot(freq,pow2db(dens_ang(:,3)),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
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


%% Spettro Campo Magnetico
xdft=fft(mag);
xdft=xdft(1:L/2+1,:);
dens_mag=(1/(2*pi*L))*abs(xdft).^2;
dens_mag(2:end-1,:)=2*dens_mag(2:end-1,:);

figure("Name","Densità Spettrale Campo Magnetico")
subplot(3,1,1)
plot(freq,pow2db(dens_mag(:,1)),LineWidth=1,Color="r");
title("Densità Spettrale Campo Magnetico")
subtitle("X''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid
subplot(3,1,2)
plot(freq,pow2db(dens_mag(:,2)),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,2])
grid
subplot(3,1,3)
plot(freq,pow2db(dens_mag(:,3)),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
% ylim([0,8])
grid


