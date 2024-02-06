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


%% Accelerazione
acc=db(inizio:fine,2:4)*9.81/-gMedio;
acc_rot=acc*gzRot;

% figure("Name","Accelerazione")
% subplot(3,1,1)
% plot(t,acc_rot(:,1),LineWidth=1,Color="r")
% title("Accelerazione")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,acc_rot(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,acc_rot(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


%% Accelerazione Media
acc_mean=movmean(acc*gzRot,[40,0]);

figure("Name","Accelerazione Media")
subplot(3,1,1)
plot(t,acc_mean(:,1),LineWidth=1,Color="r")
title("Accelerazione Media")
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


%% Velocità Angolare
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


%% Campo Magnetico
% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);

% figure("Name","Campo Magnetico")
% subplot(3,1,1)
% plot(t,mag(:,1),LineWidth=1,Color="r")
% title("Campo Magnetico")
% subtitle("X")
% xlabel("t(s)")
% ylabel("µT")
% grid
% subplot(3,1,2)
% plot(t,mag(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("µT")
% grid
% subplot(3,1,3)
% plot(t,mag(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("µT")
% grid


%% ahrs
ahrsFuse=ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');

% ahrsFuse=complementaryFilter('SampleRate',sr,'OrientationFormat','quaternion', ...
%     'ReferenceFrame','NED');

[ahrs_orient,ahrs_vang]=ahrsFuse(acc,vang,mag);


%% G
% g=rotateframe(ahrs_orient,[0,0,9.81])*gzRot;
 
% figure("Name","Gravità")
% subplot(3,1,1)
% plot(t,g(:,1),LineWidth=1,Color="r")
% title("Gravità")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,2)
% plot(t,g(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% subplot(3,1,3)
% plot(t,g(:,3),LineWidth=1,Color="b")
% subtitle("Z")
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


acc_noG=(acc-rotateframe(ahrs_orient,[0,0,9.81]))*gzRot;

figure("Name","Accelerazione Tolta la Gravità")
subplot(3,1,1)
plot(t,acc_noG(:,1),LineWidth=1,Color="r")
title("Accelerazione Tolta la Gravità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_noG(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_noG(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% ahrs Accelerazione
acc_noG_mean=movmean(acc_noG,[40,0]);

figure("Name","Media Accelerazione Tolta la Gravità")
subplot(3,1,1)
plot(t,acc_noG_mean(:,1),LineWidth=1,Color="r")
title("Media Accelerazione Tolta la Gravità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_noG_mean(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_noG_mean(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Norma
norma_acc=zeros(length(acc_rot),1);
norma_acc_media=zeros(length(acc_mean),1);
norma_acc_noG=zeros(length(acc_noG),1);
norma_acc_noG_mean=zeros(length(acc_noG_mean),1);

for i=1:length(acc_noG_mean)
    norma_acc(i)=norm(acc_rot(i,:));
    norma_acc_media(i)=norm(acc_mean(i,:));
    norma_acc_noG(i)=norm(acc_noG(i,:));
    norma_acc_noG_mean(i)=norm(acc_noG_mean(i,:));
end

% figure(Name="Norma Accelerazione")
% plot(t,norma_acc,LineWidth=1,Color="r");
% title("Norma Accelerazione")
% % ylim([0,1])
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


figure(Name="Norma Accelerazione Media")
plot(t,norma_acc_media,LineWidth=1,Color="r");
title("Norma Accelerazione Media")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


% figure(Name="Norma Accelerazione Tolta la Gravità")
% plot(t,norma_acc_noG,LineWidth=1,Color="r");
% title("Norma Accelerazione Tolta la Gravità")
% % ylim([0,1])
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


figure(Name="Norma Accelerazione Media Tolta la Gravità")
plot(t,norma_acc_noG_mean,LineWidth=1,Color="r");
title("Norma Accelerazione Media Tolta la Gravità")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Norma XY
norma_acc_XY=zeros(length(acc_rot),1);
norma_acc_media_XY=zeros(length(acc_mean),1);
norma_acc_noG_XY=zeros(length(acc_noG),1);
norma_acc_noG_mean_XY=zeros(length(acc_noG_mean),1);

for i=1:length(acc_noG_mean)
    norma_acc_XY(i)=norm(acc_rot(i,1:2));
    norma_acc_media_XY(i)=norm(acc_mean(i,1:2));
    norma_acc_noG_XY(i)=norm(acc_noG(i,1:2));
    norma_acc_noG_mean_XY(i)=norm(acc_noG_mean(i,1:2));
end

% figure(Name="Norma Accelerazione XY")
% plot(t,norma_acc_XY,LineWidth=1,Color="r");
% title("Norma Accelerazione XY")
% % ylim([0,1])
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


figure(Name="Norma Accelerazione Media XY")
plot(t,norma_acc_media_XY,LineWidth=1,Color="r");
title("Norma Accelerazione Media XY")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid


% figure(Name="Norma Accelerazione Tolta la Gravità XY")
% plot(t,norma_acc_noG_XY,LineWidth=1,Color="r");
% title("Norma Accelerazione Tolta la Gravità XY")
% % ylim([0,1])
% xlabel("t(s)")
% ylabel("m/s^2")
% grid


figure(Name="Norma Accelerazione Media Tolta la Gravità XY")
plot(t,norma_acc_noG_mean_XY,LineWidth=1,Color="r");
title("Norma Accelerazione Media Tolta la Gravità XY")
% ylim([0,1])
xlabel("t(s)")
ylabel("m/s^2")
grid



