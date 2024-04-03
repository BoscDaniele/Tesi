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
rilievo=4;


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
acc_rotf=lowpass(acc,0.5,sr)*gzRot;

figure("Name","Accelerazione")
subplot(3,1,1)
plot(t,acc_rotf(:,1),LineWidth=1,Color="r")
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

% figure("Name","Velocità Angolare")
% subplot(3,1,1)
% plot(t,vang(:,1),LineWidth=1,Color="r")
% title("Velocità Angolare")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,2)
% plot(t,vang(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,3)
% plot(t,vang(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad/s")
% grid


% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);


%% Accelerazione
acc_mean=movmean(acc_rot,[40,0]);
acc_median=movmedian(acc_rot,[40,0]);
acc_var=movvar(acc_rot,[40,0]);
acc_std=movstd(acc_rot,[40,0]);
acc_max=movmax(acc_rot,[40,0]);

acc_picco=peak2peak(acc_rot(:,1));


figure("Name","Media Accelerazione")
subplot(3,1,1)
plot(t,acc_mean(:,1),LineWidth=1,Color="r")
title("media Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_mean(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_mean(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


% figure("Name","Mediana Accelerazione")
% subplot(3,1,1)
% plot(t,acc_median(:,1),LineWidth=1,Color="r")
% title("mediana Accelerazione")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_median(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_median(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


% figure("Name","Varianza Accelerazione")
% subplot(3,1,1)
% plot(t,acc_var(:,1),LineWidth=1,Color="r")
% title("Varianza Accelerazione")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_var(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


% figure("Name","Deviazione Standard Accelerazione")
% subplot(3,1,1)
% plot(t,acc_std(:,1),LineWidth=1,Color="r")
% title("Deviazione Standard Accelerazione")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_std(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_std(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


% figure("Name","Max Accelerazione")
% subplot(3,1,1)
% plot(t,acc_max(:,1),LineWidth=1,Color="r")
% title("Max Accelerazione")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_max(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_max(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


% figure("Name","Accelerazione Picco-Picco")
% plot(t,acc_picco,LineWidth=1,Color="r")
% title("Accelerazione Picco-Picco")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


%% Velocità
vel=cumsum(acc*gzRot)*0.04;
vel_media=movmean(vel,[40,0]);
vel_var=movvar(vel,[40,0]);


% figure("Name","Media Velocità")
% subplot(3,1,1)
% plot(t,vel_media(:,1),LineWidth=1,Color="r")
% title("media Velocità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,2)
% plot(t,vel_media(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,3)
% plot(t,vel_media(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s")
% grid


% figure("Name","Varianza Velocità")
% subplot(3,1,1)
% plot(t,vel_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,2)
% plot(t,vel_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,3)
% plot(t,vel_var(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s")
% grid


%% Velocità Angolare
vang_mean=movmean(vang,[40,0]);
vang_var=movvar(vang,[40,0]);


% figure("Name","Media Velocità Angolare")
% subplot(3,1,1)
% plot(t,vang_mean(:,1),LineWidth=1,Color="r")
% title("Media Velocità Angolare")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,2)
% plot(t,vang_mean(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,3)
% plot(t,vang_mean(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad/s")
% grid


% figure("Name","Varianza Velocità Angolare")
% subplot(3,1,1)
% plot(t,vang_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,2)
% plot(t,vang_var(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,3)
% plot(t,vang_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad/s")
% grid


%% Integrale Velocità Angolare
ang=cumsum(vang);
ang_media=movmean(ang,[40,0]);
ang_var=movvar(ang,[40,0]);


% figure("Name","Integrale Velocità Angolare Medio")
% subplot(3,1,1)
% plot(t,ang_media(:,1),LineWidth=1,Color="r")
% title("Integrale Velocità Angolare Medio")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,2)
% plot(t,ang_media(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,3)
% plot(t,ang_media(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad")
% grid


% figure("Name","Varianza Integrale Velocità Angolare")
% subplot(3,1,1)
% plot(t,ang_var(:,1),LineWidth=1,Color="r")
% title("Varianza Integrale Velocità Angolare")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,2)
% plot(t,ang_var(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,3)
% plot(t,ang_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad")
% grid


% %% ahrs
% ahrsFuse=ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
%     'ReferenceFrame','NED');
% 
% [ahrs_orient,ahrs_vang]=ahrsFuse(acc,vang,mag);
% 
% acc_noG=(acc-rotateframe(ahrs_orient,[0,0,9.81]))*gzRot;
% 
% 
% %% ahrs Accelerazione
% acc_noG_mean=movmean(acc_noG,[40,0]);
% acc_noG_var=movvar(acc_noG,[40,0]);
% 
% figure("Name","Media Accelerazione Tolta la Gravità")
% subplot(3,1,1)
% plot(t,acc_noG_mean(:,1),LineWidth=1,Color="r")
% title("Media Accelerazione Tolta la Gravità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_noG_mean(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_noG_mean(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% 
% 
% figure("Name","Varianza Accelerazione Tolta la Gravità")
% subplot(3,1,1)
% plot(t,acc_noG_var(:,1),LineWidth=1,Color="r")
% title("Varianza Accelerazione Tolta la Gravità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_noG_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_noG_var(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% 
% 
% %% ahrs Velocità
% vel_noG=cumsum(acc_noG)*0.04;
% vel_noG_media=movmean(vel_noG,[40,0]);
% vel_noG_var=movvar(vel_noG,[40,0]);
% 
% 
% figure("Name","Media Velocità Tolta la Gravità")
% subplot(3,1,1)
% plot(t,vel_noG_media(:,1),LineWidth=1,Color="r")
% title("Media Velocità Tolta la Gravità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,2)
% plot(t,vel_noG_media(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,3)
% plot(t,vel_noG_media(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% 
% 
% figure("Name","Varianza Velocità Tolta la Gravità")
% subplot(3,1,1)
% plot(t,vel_noG_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Tolta la Gravità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,2)
% plot(t,vel_noG_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% subplot(3,1,3)
% plot(t,vel_noG_var(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s")
% grid
% 
% 
% %% ahrs Velocità Angolare
% ahrs_vang_media=movmean(ahrs_vang,[40,0]);
% ahrs_vang_var=movvar(ahrs_vang,[40,0]);
% 
% 
% figure("Name","Media Velocità Angolare AHRS")
% subplot(3,1,1)
% plot(t,ahrs_vang_media(:,1),LineWidth=1,Color="r")
% title("Media Velocità Angolare AHRS")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,2)
% plot(t,ahrs_vang_media(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,3)
% plot(t,ahrs_vang_media(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% 
% 
% figure("Name","Varianza Velocità Angolare AHRS")
% subplot(3,1,1)
% plot(t,ahrs_vang_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare AHRS")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,2)
% plot(t,ahrs_vang_var(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% subplot(3,1,3)
% plot(t,ahrs_vang_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad/s")
% grid
% 
% 
% %% ahrs Integrale Velocità Angolare
% ahrs_ang=cumsum(ahrs_vang)*0.04;
% ahrs_ang_media=movmean(ahrs_ang,[40,0]);
% ahrs_ang_var=movvar(ahrs_ang,[40,0]);
% 
% 
% figure("Name","Integrale Velocità Angolare Medio AHRS")
% subplot(3,1,1)
% plot(t,ahrs_ang_media(:,1),LineWidth=1,Color="r")
% title("Integrale Velocità Angolare Medio AHRS")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,2)
% plot(t,ahrs_ang_media(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,3)
% plot(t,ahrs_ang_media(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad")
% grid
% 
% 
% figure("Name","Varianza Integrale Velocità Angolare AHRS")
% subplot(3,1,1)
% plot(t,ahrs_ang_var(:,1),LineWidth=1,Color="r")
% title("Varianza Integrale Velocità Angolare AHRS")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,2)
% plot(t,ahrs_ang_var(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("rad")
% grid
% subplot(3,1,3)
% plot(t,ahrs_ang_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("rad")
% grid
% 
% 
% 
% 
% 



