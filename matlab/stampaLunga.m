clear all
close all
clc

pathP=".\dati\lunga_piano\";
pathF=".\dati\lunga_forte\";

rilievoP=3;
rilievoF=2;

[gzRotP,gMedioP] = GZRot(pathP);
[gzRotF,gMedioF] = GZRot(pathF);

sr = 25; %sample rate

%% Dati Piano
dbP=importdata(pathP + "BlueCoin_Log_N00"+rilievoP+".csv").data;

tP_rotta=dbP(:,1)*1e-3;
tP_rotta=tP_rotta-tP_rotta(1);

accP_rotta=dbP(:,2:4)*gzRotP*9.81/-gMedioP;
vangP_rotta=(dbP(:,5:7)*1e-3);
magP_rotta=([dbP(:,8),-dbP(:,9),dbP(:,10)]*1e-1);

[tP,accP,vangP,magP]=AggiustaFrequenza(tP_rotta,accP_rotta,vangP_rotta,magP_rotta);


%% Dati Forte
dbF=importdata(pathF + "BlueCoin_Log_N00"+rilievoF+".csv").data;

tF_rotta=dbF(:,1)*1e-3;
tF_rotta=tF_rotta-tF_rotta(1);

accF_rotta=dbF(:,2:4)*gzRotF*9.81/-gMedioF;
vangF_rotta=(dbF(:,5:7)*1e-3);
magF_rotta=([dbF(:,8),-dbF(:,9),dbF(:,10)]*1e-1);

[tF,accF,vangF,magF]=AggiustaFrequenza(tF_rotta,accF_rotta,vangF_rotta,magF_rotta);


%%
accP_media=movmean(accP,40);
accP_mediaPost=movmean(accP,[40,0]);
accP_lowpass=lowpass(accP,0.5,25);

accF_media=movmean(accF,40);
accF_mediaPost=movmean(accF,[40,0]);
accF_lowpass=lowpass(accF,0.5,25);



% %% Accelerazione XY Forte/Piano
% limX=[floor(min(min(accP(:,1)),min(accF(:,1)))),ceil(max(max(accP(:,1)),max(accF(:,1))))];
% limY=[floor(min(min(accP(:,2)),min(accF(:,2)))),ceil(max(max(accP(:,2)),max(accF(:,2))))];
% 
% figure
% subplot(2,2,1)
% plot(tP,accP(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid
% subplot(2,2,2)
% plot(tF,accF(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid


% %% LowPass Accelerazione
% limX=[floor(min(min(accP_lowpass(:,1)),min(accF_lowpass(:,1)))),ceil(max(max(accP_lowpass(:,1)),max(accF_lowpass(:,1))))];
% 
% figure
% subplot(2,1,1)
% plot(tP,accP_lowpass(:,1),LineWidth=1,Color="r")
% title("Passa Basso Accelerazione Piano")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tP,zeros(length(tP)),LineWidth=.1,Color="black")
% 
% subplot(2,1,2)
% plot(tF,accF_lowpass(:,1),LineWidth=1,Color="r")
% title("Passa Basso Accelerazione Forte")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tF,zeros(length(tF)),LineWidth=.1,Color="black")


% %% Accelerazione X - Media
% figure
% plot(tF,accF(:,1),LineWidth=1)
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% hold on
% plot(tF,accF_media(:,1),LineWidth=1)
% plot(tF,zeros(length(tF),1),LineWidth=.1,Color="black")
% legend("Accelerazione X", "Media Accelerazione X")


% %% Media vs Media Post
% figure
% plot(tF,accF_media(:,1),LineWidth=1)
% xlabel("t(s)")
% ylabel("m/s^2")
% grid
% hold on
% plot(tF,accF_mediaPost(:,1),LineWidth=1)
% legend("Media Accelerazione","Media Accelerazione Postuma")


% %% Trasformata e Spettro Accelerazione
% 
% % piano
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% YP=fft(accP);
% P2P=abs(YP/LP);
% trasform_accP=P2P(1:(LP/2+1),:);
% trasform_accP(2:end-1,:)=2*trasform_accP(2:end-1,:);
% 
% xdftP=YP(1:LP/2+1,:);
% spettro_accP=(1/(sr*LP))*abs(xdftP).^2;
% spettro_accP(2:end-1,:)=2*spettro_accP(2:end-1,:);
% 
% % forte
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% YF=fft(accF);
% P2F=abs(YF/LF);
% trasform_accF=P2F(1:(LF/2+1),:);
% trasform_accF(2:end-1,:)=2*trasform_accF(2:end-1,:);
% 
% xdftF=YF(1:LF/2+1,:);
% spettro_accF=(1/(sr*LF))*abs(xdftF).^2;
% spettro_accF(2:end-1,:)=2*spettro_accF(2:end-1,:);
% 
% % Trasformata Accelerazione X
% figure
% plot(fP,trasform_accP(:,1),LineWidth=1,Color="r")
% title("Trasformata Accelerazione X")
% xlabel("f(Hz)")
% ylabel("X''(f)")
% grid
% hold on
% plot(fF,trasform_accF(:,1),LineWidth=1,Color="b")
% legend("Trasformata Acc Piano","Trasformata Acc Forte")
% 
% % Trasformata Accelerazione Y
% figure
% plot(fP,trasform_accP(:,2),LineWidth=1,Color="r")
% title("Trasformata Accelerazione Y")
% xlabel("f(Hz)")
% ylabel("X''(f)")
% grid
% hold on
% plot(fF,trasform_accF(:,2),LineWidth=1,Color="b")
% legend("Trasformata Acc Piano","Trasformata Acc Forte")
% 
% 
% % Spettro Accelerazione X
% figure
% plot(fP,spettro_accP(:,1),LineWidth=1,Color="r")
% title("Spettro Accelerazione X")
% xlabel("f(Hz)")
% ylabel("X''(f)")
% grid
% hold on
% plot(fF,spettro_accF(:,1),LineWidth=1,Color="b")
% legend("Spettro Acc Piano","Spettro Acc Forte")
% 
% % Spettro Accelerazione Y
% figure
% plot(fP,spettro_accP(:,2),LineWidth=1,Color="r")
% title("Spettro Accelerazione Y")
% xlabel("f(Hz)")
% ylabel("X''(f)")
% grid
% hold on
% plot(fF,spettro_accF(:,2),LineWidth=1,Color="b")
% legend("Spettro Acc Piano","Spettro Acc Forte")


% %% Media Accelerazione
% limX=[floor(min(min(accP_media(:,1)),min(accF_media(:,1)))),ceil(max(max(accP_media(:,1)),max(accF_media(:,1))))];
% 
% figure
% subplot(2,1,1)
% plot(tP,accP_media(:,1),LineWidth=1,Color="r")
% title("Media Accelerazione Piano")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tP,zeros(length(tP)),LineWidth=.1,Color="black")
% 
% subplot(2,1,2)
% plot(tF,accF_media(:,1),LineWidth=1,Color="r")
% title("Media Accelerazione Forte")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tF,zeros(length(tF)),LineWidth=.1,Color="black")


% %% Media Posticipata Accelerazione
% limX=[floor(min(min(accP_mediaPost(:,1)),min(accF_mediaPost(:,1)))),ceil(max(max(accP_mediaPost(:,1)),max(accF_mediaPost(:,1))))];
% 
% figure
% subplot(2,1,1)
% plot(tP,accP_mediaPost(:,1),LineWidth=1,Color="r")
% title("Media Posticipata Accelerazione Piano")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tP,zeros(length(tP)),LineWidth=.1,Color="black")
% 
% subplot(2,1,2)
% plot(tF,accF_mediaPost(:,1),LineWidth=1,Color="r")
% title("Media Posticipata Accelerazione Forte")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% hold on
% plot(tF,zeros(length(tF)),LineWidth=.1,Color="black")


% %% Varianza Accelerazione
% accP_var=movvar(accP,40);
% accF_var=movvar(accF,40);
% 
% limX=[floor(min(min(accP_var(:,1)),min(accF_var(:,1)))),ceil(max(max(accP_var(:,1)),max(accF_var(:,1))))];
% limY=[floor(min(min(accP_var(:,2)),min(accF_var(:,2)))),ceil(max(max(accP_var(:,2)),max(accF_var(:,2))))];
% 
% figure
% subplot(2,2,1)
% plot(tP,accP_var(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_var(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_var(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid


% %% Deviazione Standard Accelerazione
% accP_std=movstd(accP,40);
% accF_std=movstd(accF,40);
% 
% limX=[floor(min(min(accP_std(:,1)),min(accF_std(:,1)))),ceil(max(max(accP_std(:,1)),max(accF_std(:,1))))];
% limY=[floor(min(min(accP_std(:,2)),min(accF_std(:,2)))),ceil(max(max(accP_std(:,2)),max(accF_std(:,2))))];
% 
% figure
% subplot(2,2,1)
% plot(tP,accP_std(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_std(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_std(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte")
% subtitle("X")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_std(:,2),LineWidth=1,Color="g")
% subtitle("Y")
% xlabel("t(s)")
% ylabel("m/s^2")
% ylim(limY)
% grid


% %% Velocità Angolare
% 
% limX=[floor(min(min(vangP(:,1)),min(vangF(:,1)))),ceil(max(max(vangP(:,1)),max(vangF(:,1))))];
% limY=[floor(min(min(vangP(:,2)),min(vangF(:,2)))),ceil(max(max(vangP(:,2)),max(vangF(:,2))))];
% limZ=[floor(min(min(vangP(:,3)),min(vangF(:,3)))),ceil(max(max(vangP(:,3)),max(vangF(:,3))))];
% 
% figure
% subplot(3,2,1)
% plot(tP,vangP(:,1),LineWidth=1,Color="r")
% title("Velocità Angolare Piano")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limX)
% grid
% subplot(3,2,3)
% plot(tP,vangP(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limY)
% grid
% subplot(3,2,5)
% plot(tP,vangP(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limZ)
% grid
% 
% subplot(3,2,2)
% plot(tF,vangF(:,1),LineWidth=1,Color="r")
% title("Velocità Angolare Forte")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limX)
% grid
% subplot(3,2,4)
% plot(tF,vangF(:,2),LineWidth=1,Color="g")
% subtitle("Pitch")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limY)
% grid
% subplot(3,2,6)
% plot(tF,vangF(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limZ)
% grid


% %% Varianza Velocità Angolare
% vangP_var=movvar(vangP,40);
% vangF_var=movvar(vangF,40);
% 
% 
% limX=[floor(min(min(vangP_var(:,1)),min(vangF_var(:,1)))),ceil(max(max(vangP_var(:,1)),max(vangF_var(:,1))))];
% limZ=[floor(min(min(vangP_var(:,3)),min(vangF_var(:,3)))),ceil(max(max(vangP_var(:,3)),max(vangF_var(:,3))))];
% 
% figure
% subplot(2,2,1)
% plot(tP,vangP_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare Piano")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,vangP_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limZ)
% grid
% 
% subplot(2,2,2)
% plot(tF,vangF_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare Forte")
% subtitle("Roll")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,vangF_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw")
% xlabel("t(s)")
% ylabel("deg/s")
% ylim(limZ)
% grid


% %% Velocità
% velP=cumsum(accP(:,1));
% velF=cumsum(accF(:,1));
% 
% lim=[floor(min(min(velP),min(velF))),ceil(max(max(velP),max(velF)))];
% 
% figure
% subplot(2,1,1)
% plot(tP,velP,LineWidth=1)
% title("Velocità")
% subtitle("Piano")
% xlabel("t(s)")
% ylabel("m/s")
% ylim(lim)
% grid
% hold on
% plot(tP,zeros(length(velP)),LineWidth=.1,Color="black")
% subplot(2,1,2)
% plot(tF,velF,LineWidth=1)
% subtitle("Forte")
% xlabel("t(s)")
% ylabel("m/s")
% ylim(lim)
% grid
% hold on
% plot(tF,zeros(length(velF)),LineWidth=.1,Color="black")









