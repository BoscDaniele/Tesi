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


%% Trasformata Accelerazione
inizio_f=1;

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


%% Trasformata Velocità
Y=fft(vel);
P2=abs(Y/L);
trasform_vel=P2(1:(L/2+1),:);
trasform_vel(2:end-1,:)=2*trasform_vel(2:end-1,:);

figure("Name","Trasformata Velocità")
subplot(3,1,1)
plot(f(inizio_f:end),trasform_vel(inizio_f:end,1),LineWidth=1,Color="r");
title("Trasformata Velocità")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,1])
grid
subplot(3,1,2)
plot(f(inizio_f:end),trasform_vel(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,1])
grid
subplot(3,1,3)
plot(f(inizio_f:end),trasform_vel(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,1])
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


%% Trasformata Integrale Velocità Angolare
Y=fft(ang);
P2=abs(Y/L);
trasform_int_vang=P2(1:(L/2+1),:);
trasform_int_vang(2:end-1,:)=2*trasform_int_vang(2:end-1,:);

figure("Name","Trasformata Integrale Velocità Angolare")
subplot(3,1,1)
plot(f(inizio_f:end),trasform_int_vang(inizio_f:end,1),LineWidth=1,Color="r");
title("Trasformata Integrale Velocità Angolare")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,1])
grid
subplot(3,1,2)
plot(f(inizio_f:end),trasform_int_vang(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,1])
grid
subplot(3,1,3)
plot(f(inizio_f:end),trasform_int_vang(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,1])
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


%% Trasformata Campo Magnetico
Y=fft(mag);
P2=abs(Y/L);
trasform_mag=P2(1:(L/2+1),:);
trasform_mag(2:end-1,:)=2*trasform_mag(2:end-1,:);

figure("Name","Trasformata Campo Magnetico")
subplot(3,1,1)
plot(f(inizio_f:end),trasform_mag(inizio_f:end,1),LineWidth=1,Color="r");
title("Trasformata Campo Magnetico")
subtitle("X''(f)")
xlabel("f(Hz)")
ylabel("X''(f)")
% ylim([0,1])
grid
subplot(3,1,2)
plot(f(inizio_f:end),trasform_mag(inizio_f:end,2),LineWidth=1,Color="g");
subtitle("Y''(f)")
xlabel("f(Hz)")
ylabel("Y''(f)")
% ylim([0,1])
grid
subplot(3,1,3)
plot(f(inizio_f:end),trasform_mag(inizio_f:end,3),LineWidth=1,Color="b");
subtitle("Z''(f)")
xlabel("f(Hz)")
ylabel("Z''(f)")
% ylim([0,1])
grid
