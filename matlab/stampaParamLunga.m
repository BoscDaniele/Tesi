clear all
close all
clc

path=".\dati\lunga_forte\";

rilievo=2;

[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*gzRot*9.81/-gMedio;
vang_rotta=(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1);

[t,acc,vang,mag]=AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta);


% %% Accelerazione
% stampa(t,acc,"Accelerazione",['X','Y'],'t(s)','m/s^2')


% %% LowPass Accelerazione
% acc_low=lowpass(acc,0.5,25);
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_low(:,i))),ceil(max(acc_low(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2)+0.5;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_low(:,1:2),"AccLowPass",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Media Accelerazione
% acc_media=movmean(acc,40);
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_media(:,i))),ceil(max(acc_media(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2);
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_media(:,1:2),"AccMedia",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Valore Medio Rettificato
% acc_arv=movmean(abs(acc),40);
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_arv(:,i))),ceil(max(acc_arv(:,i)))];
% end
% limY(1,2)=limY(1,2)+0.5;
% limY(2,2)=limY(2,2);
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_arv(:,1:2),"AccMediaRett",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Varianza Accelerazione
% acc_var=movvar(acc,40);
% stampa(t,acc_var,"AccVar",['X','Y'],'t(s)','m/s^2')


%% Deviazione Standard Accelerazione
acc_std=movstd(acc,40);

limY=ones(3,2);
for i=1:3
    limY(i,:)=[floor(min(acc_std(:,i))),ceil(max(acc_std(:,i)))];
end
limY(1,2)=limY(1,2)+1;
limY(2,2)=limY(2,2);

textLimY=ones(2,5);
textLimY(1,:)=(limY(1,2)-0.5)*ones(5,1);
textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);

stampa_gen(t,acc_std(:,1:2),"AccStd",['X','Y'],'t(s)','m/s^2',limY,textLimY)


%% Scarto Quadratico Medio Accelerazione
acc_rms=movrms(acc,40);

limY=ones(3,2);
for i=1:3
    limY(i,:)=[floor(min(acc_rms(:,i))),ceil(max(acc_rms(:,i)))];
end
limY(1,2)=limY(1,2)+0.5;
limY(2,2)=limY(2,2);

textLimY=ones(2,5);
textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);

stampa_gen(t,acc_rms(:,1:2),"AccRms",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Kurtosi Accelerazione
% accP_krt=movkurt(accP,40);
% accF_krt=movkurt(accF,40);
% 
% limX=[floor(min(min(accP_krt(:,1)),min(accF_krt(:,1)))),ceil(max(max(accP_krt(:,1)),max(accF_krt(:,1))))];
% limY=[floor(min(min(accP_krt(:,2)),min(accF_krt(:,2)))),ceil(max(max(accP_krt(:,2)),max(accF_krt(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_krt(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_krt(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_krt(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_krt(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_krtAccXY.png")
% end


% %% Skewness Accelerazione
% accP_skw=movskw(accP,40);
% accF_skw=movskw(accF,40);
% 
% limX=[floor(min(min(accP_skw(:,1)),min(accF_skw(:,1)))),ceil(max(max(accP_skw(:,1)),max(accF_skw(:,1))))];
% limY=[floor(min(min(accP_skw(:,2)),min(accF_skw(:,2)))),ceil(max(max(accP_skw(:,2)),max(accF_skw(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_skw(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_skw(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_skw(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_skw(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_skwAccXY.png")
% end


% %% Shape Factor Accelerazione 
% accP_rms=movrms(accP,40);
% accF_rms=movrms(accF,40);
% 
% accP_arv=movmean(abs(accP),[40,0]);
% accF_arv=movmean(abs(accF),[40,0]);
% 
% accP_shf=accP_rms./accP_arv;
% accF_shf=accF_rms./accF_arv;
% 
% limX=[floor(min(min(accP_shf(:,1)),min(accF_shf(:,1)))),ceil(max(max(accP_shf(:,1)),max(accF_shf(:,1))))];
% limY=[floor(min(min(accP_shf(:,2)),min(accF_shf(:,2)))),ceil(max(max(accP_shf(:,2)),max(accF_shf(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_shf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_shf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_shf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_shf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_shfAccXY.png")
% end


% %% Crest Factor Accelerazione
% accP_max=movmax(accP,[40,0]);
% accF_max=movmax(accF,[40,0]);
% 
% accP_rms=movrms(accP,40);
% accF_rms=movrms(accF,40);
% 
% accP_crf=accP_max./accP_rms;
% accF_crf=accF_max./accF_rms;
% 
% limX=[floor(min(min(accP_crf(:,1)),min(accF_crf(:,1)))),ceil(max(max(accP_crf(:,1)),max(accF_crf(:,1))))];
% limY=[floor(min(min(accP_crf(:,2)),min(accF_crf(:,2)))),ceil(max(max(accP_crf(:,2)),max(accF_crf(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_crf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_crf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_crf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_crf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_crfAccXY.png")
% end


% %% Impulse Factor Accelerazione
% accP_max=movmax(accP,[40,0]);
% accF_max=movmax(accF,[40,0]);
% 
% accP_arv=movmean(abs(accP),[40,0]);
% accF_arv=movmean(abs(accF),[40,0]);
% 
% accP_impf=accP_max./accP_arv;
% accF_impf=accF_max./accF_arv;
% 
% limX=[floor(min(min(accP_impf(:,1)),min(accF_impf(:,1)))),ceil(max(max(accP_impf(:,1)),max(accF_impf(:,1))))];
% limY=[floor(min(min(accP_impf(:,2)),min(accF_impf(:,2)))),ceil(max(max(accP_impf(:,2)),max(accF_impf(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_impf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_impf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_impf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_impf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_impfAccXY.png")
% end


% %% Margin Factor Accelerazione
% accP_max=movmax(accP,[40,0]);
% accF_max=movmax(accF,[40,0]);
% 
% accP_arvq=movmean(sqrt(abs(accP)),[40,0]).^2;
% accF_arvq=movmean(sqrt(abs(accF)),[40,0]).^2;
% 
% accP_mrgf=accP_max./accP_arvq;
% accF_mrgf=accF_max./accF_arvq;
% 
% limX=[floor(min(min(accP_mrgf(:,1)),min(accF_mrgf(:,1)))),ceil(max(max(accP_mrgf(:,1)),max(accF_mrgf(:,1))))];
% limY=[floor(min(min(accP_mrgf(:,2)),min(accF_mrgf(:,2)))),ceil(max(max(accP_mrgf(:,2)),max(accF_mrgf(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_mrgf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_mrgf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_mrgf(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_mrgf(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_mrgfAccXY.png")
% end


% %% Max Accelerazione
% n_max=[10,0];
% 
% accP_max=movmax(accP,n_max);
% accF_max=movmax(accF,n_max);
% 
% limX=[floor(min(min(accP_max(:,1)),min(accF_max(:,1)))),ceil(max(max(accP_max(:,1)),max(accF_max(:,1))))];
% limY=[floor(min(min(accP_max(:,2)),min(accF_max(:,2)))),ceil(max(max(accP_max(:,2)),max(accF_max(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_max(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_max(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_max(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_max(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_maxAccXY.png")
% end


% %% Min Accelerazione
% n_min=[10,0];
% 
% accP_min=movmin(accP,n_min);
% accF_min=movmin(accF,n_min);
% 
% limX=[floor(min(min(accP_min(:,1)),min(accF_min(:,1)))),ceil(max(max(accP_min(:,1)),max(accF_min(:,1))))];
% limY=[floor(min(min(accP_min(:,2)),min(accF_min(:,2)))),ceil(max(max(accP_min(:,2)),max(accF_min(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_min(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_min(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_min(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_min(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_minAccXY.png")
% end


% %% Peak Accelerazione
% n_Peak=10;
% 
% accP_max=movmax(accP,n_Peak);
% accF_max=movmax(accF,n_Peak);
% 
% accP_min=movmin(accP,n_Peak);
% accF_min=movmin(accF,n_Peak);
% 
% accP_peak=accP_max-accP_min;
% accF_peak=accF_max-accF_min;
% 
% limX=[floor(min(min(accP_peak(:,1)),min(accF_peak(:,1)))),ceil(max(max(accP_peak(:,1)),max(accF_peak(:,1))))];
% limY=[floor(min(min(accP_peak(:,2)),min(accF_peak(:,2)))),ceil(max(max(accP_peak(:,2)),max(accF_peak(:,2))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,accP_peak(:,1),LineWidth=1,Color="r")
% title("Accelerazione Piano",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,accP_peak(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% subplot(2,2,2)
% plot(tF,accF_peak(:,1),LineWidth=1,Color="r")
% title("Accelerazione Forte",FontName=font)
% subtitle("X",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,accF_peak(:,2),LineWidth=1,Color="g")
% subtitle("Y",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s^2",FontName=font)
% ylim(limY)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_peak2peakAccXY.png")
% end


% %% Trasformata e Spettro Accelerazione
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% % trasformata
% YP=fft(accP);
% P2P=abs(YP/LP);
% trasform_accP=P2P(1:(LP/2+1),:);
% trasform_accP(2:end-1,:)=2*trasform_accP(2:end-1,:);
% 
% YF=fft(accF);
% P2F=abs(YF/LF);
% trasform_accF=P2F(1:(LF/2+1),:);
% trasform_accF(2:end-1,:)=2*trasform_accF(2:end-1,:);
% 
% % trasformata no media
% YP_noMedia=fft(accP-accP_media);
% P2P_noMedia=abs(YP_noMedia/LP);
% trasform_accP_noMedia=P2P_noMedia(1:(LP/2+1),:);
% trasform_accP_noMedia(2:end-1,:)=2*trasform_accP_noMedia(2:end-1,:);
% 
% YF_noMedia=fft(accF-accF_media);
% P2F_noMedia=abs(YF_noMedia/LF);
% trasform_accF_noMedia=P2F_noMedia(1:(LF/2+1),:);
% trasform_accF_noMedia(2:end-1,:)=2*trasform_accF_noMedia(2:end-1,:);
% 
% % trasformata no media/varianza
% YP_norm=fft((accP-accP_media)./movvar(accP,40));
% P2P_norm=abs(YP_norm/LP);
% trasform_accP_norm=P2P_norm(1:(LP/2+1),:);
% trasform_accP_norm(2:end-1,:)=2*trasform_accP_norm(2:end-1,:);
% 
% YF_norm=fft((accF-accF_media)./movvar(accF,40));
% P2F_norm=abs(YF_norm/LF);
% trasform_accF_norm=P2F_norm(1:(LF/2+1),:);
% trasform_accF_norm(2:end-1,:)=2*trasform_accF_norm(2:end-1,:);
% 
% 
% figure(Name="Trasformata")
% plot(fF,trasform_accF(:,1),LineWidth=1,Color="b");
% title("Trasformata X",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,trasform_accP(:,1),LineWidth=1,Color="r")
% legend("Trasformata Acc Forte","Trasformata Acc Piano",FontName=font)
% 
% 
% figure(Name="Trasformata No Media")
% plot(fF,trasform_accF_noMedia(:,1),LineWidth=1,Color="b");
% title("Trasformata No Media X",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,trasform_accP_noMedia(:,1),LineWidth=1,Color="r")
% legend("Trasformata Acc Forte","Trasformata Acc Piano",FontName=font)
% 
% 
% figure(Name="Trasformata Normalizzata")
% plot(fF,trasform_accF_norm(:,1),LineWidth=1,Color="b");
% title("Trasformata (Acc X - Media Acc X)/Varianza",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,trasform_accP_norm(:,1),LineWidth=1,Color="r")
% legend("Trasformata Acc Forte","Trasformata Acc Piano",FontName=font)


% %% Ampiezza Media
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% % trasformata
% YP=fft(accP);
% P2P=abs(YP/LP);
% trasform_accP=P2P(1:(LP/2+1),:);
% trasform_accP(2:end-1,:)=2*trasform_accP(2:end-1,:);
% 
% YF=fft(accF);
% P2F=abs(YF/LF);
% trasform_accF=P2F(1:(LF/2+1),:);
% trasform_accF(2:end-1,:)=2*trasform_accF(2:end-1,:);
% 
% average_trasform_accP=mean(trasform_accP);
% average_trasform_accF=mean(trasform_accF);
% 
% % trasformata no media
% YP_noMedia=fft(accP-accP_media);
% P2P_noMedia=abs(YP_noMedia/LP);
% trasform_accP_noMedia=P2P_noMedia(1:(LP/2+1),:);
% trasform_accP_noMedia(2:end-1,:)=2*trasform_accP_noMedia(2:end-1,:);
% 
% YF_noMedia=fft(accF-accF_media);
% P2F_noMedia=abs(YF_noMedia/LF);
% trasform_accF_noMedia=P2F_noMedia(1:(LF/2+1),:);
% trasform_accF_noMedia(2:end-1,:)=2*trasform_accF_noMedia(2:end-1,:);
% 
% average_trasform_accP_noMedia=mean(trasform_accP_noMedia);
% average_trasform_accF_noMedia=mean(trasform_accF_noMedia);
% 
% figure(Name="Trasformata")
% subplot(2,1,1)
% plot(fP,trasform_accP(:,1),LineWidth=1,Color="b");
% title("Trasformata Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,average_trasform_accP(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Ampiezza Media",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF(:,1),LineWidth=1,Color="r")
% title("Trasformata Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,average_trasform_accF(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Ampiezza Media",FontName=font)
% 
% 
% figure(Name="Trasformata No Media")
% subplot(2,1,1)
% plot(fP,trasform_accP_noMedia(:,1),LineWidth=1,Color="b");
% title("Trasformata No Media Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,average_trasform_accP_noMedia(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Ampiezza Media",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF_noMedia(:,1),LineWidth=1,Color="r")
% title("Trasformata No Media Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,average_trasform_accF_noMedia(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Ampiezza Media",FontName=font)


% %% Frequency Centroid
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% % trasformata
% YP=fft(accP);
% P2P=abs(YP/LP);
% trasform_accP=P2P(1:(LP/2+1),:);
% trasform_accP(2:end-1,:)=2*trasform_accP(2:end-1,:);
% 
% YF=fft(accF);
% P2F=abs(YF/LF);
% trasform_accF=P2F(1:(LF/2+1),:);
% trasform_accF(2:end-1,:)=2*trasform_accF(2:end-1,:);
% 
% centroid_trasform_accP=mean(trasform_accP.*fP');
% centroid_trasform_accF=mean(trasform_accF.*fF');
% 
% % trasformata no media
% YP_noMedia=fft(accP-accP_media);
% P2P_noMedia=abs(YP_noMedia/LP);
% trasform_accP_noMedia=P2P_noMedia(1:(LP/2+1),:);
% trasform_accP_noMedia(2:end-1,:)=2*trasform_accP_noMedia(2:end-1,:);
% 
% YF_noMedia=fft(accF-accF_media);
% P2F_noMedia=abs(YF_noMedia/LF);
% trasform_accF_noMedia=P2F_noMedia(1:(LF/2+1),:);
% trasform_accF_noMedia(2:end-1,:)=2*trasform_accF_noMedia(2:end-1,:);
% 
% centroid_trasform_accP_noMedia=mean(trasform_accP_noMedia.*fP');
% centroid_trasform_accF_noMedia=mean(trasform_accF_noMedia.*fF');
% 
% figure(Name="Trasformata")
% subplot(2,1,1)
% plot(fP,trasform_accP(:,1),LineWidth=1,Color="b");
% title("Trasformata Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,centroid_trasform_accP(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Centroid",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF(:,1),LineWidth=1,Color="r")
% title("Trasformata Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,centroid_trasform_accF(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Centroid",FontName=font)
% 
% 
% figure(Name="Trasformata No Media")
% subplot(2,1,1)
% plot(fP,trasform_accP_noMedia(:,1),LineWidth=1,Color="b");
% title("Trasformata No Media Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,centroid_trasform_accP_noMedia(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Centroid",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF_noMedia(:,1),LineWidth=1,Color="r")
% title("Trasformata No Media Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,centroid_trasform_accF_noMedia(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Centroid",FontName=font)


% %% Frequency Variance
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% % trasformata
% YP=fft(accP);
% P2P=abs(YP/LP);
% trasform_accP=P2P(1:(LP/2+1),:);
% trasform_accP(2:end-1,:)=2*trasform_accP(2:end-1,:);
% 
% YF=fft(accF);
% P2F=abs(YF/LF);
% trasform_accF=P2F(1:(LF/2+1),:);
% trasform_accF(2:end-1,:)=2*trasform_accF(2:end-1,:);
% 
% centroidP=mean(trasform_accP.*fP');
% centroidF=mean(trasform_accF.*fF');
% 
% 
% var_trasform_accP=sum((fP'-centroidP).*trasform_accP(:,1))/sum(trasform_accP(:,1));
% var_trasform_accF=sum((fF'-centroidF).*trasform_accF(:,1))/sum(trasform_accF(:,1));
% 
% % trasformata no media
% YP_noMedia=fft(accP-accP_media);
% P2P_noMedia=abs(YP_noMedia/LP);
% trasform_accP_noMedia=P2P_noMedia(1:(LP/2+1),:);
% trasform_accP_noMedia(2:end-1,:)=2*trasform_accP_noMedia(2:end-1,:);
% 
% YF_noMedia=fft(accF-accF_media);
% P2F_noMedia=abs(YF_noMedia/LF);
% trasform_accF_noMedia=P2F_noMedia(1:(LF/2+1),:);
% trasform_accF_noMedia(2:end-1,:)=2*trasform_accF_noMedia(2:end-1,:);
% 
% centroidP_noMedia=mean(trasform_accP_noMedia.*fP');
% centroidF_noMedia=mean(trasform_accF_noMedia.*fF');
% 
% 
% var_trasform_accP_noMedia=sum((fP'-centroidP_noMedia).*trasform_accP_noMedia(:,1))/sum(trasform_accP_noMedia(:,1));
% var_trasform_accF_noMedia=sum((fF'-centroidF_noMedia).*trasform_accF_noMedia(:,1))/sum(trasform_accF_noMedia(:,1));
% 
% 
% figure(Name="Trasformata")
% subplot(2,1,1)
% plot(fP,trasform_accP(:,1),LineWidth=1,Color="b");
% title("Trasformata Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,var_trasform_accP(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Variance",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF(:,1),LineWidth=1,Color="r")
% title("Trasformata Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,var_trasform_accF(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Variance",FontName=font)
% 
% 
% figure(Name="Trasformata No Media")
% subplot(2,1,1)
% plot(fP,trasform_accP_noMedia(:,1),LineWidth=1,Color="b");
% title("Trasformata No Media Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,var_trasform_accP_noMedia(1)*ones(length(trasform_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Variance",FontName=font)
% subplot(2,1,2)
% plot(fF,trasform_accF_noMedia(:,1),LineWidth=1,Color="r")
% title("Trasformata No Media Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,var_trasform_accF_noMedia(1)*ones(length(trasform_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Frequency Variance",FontName=font)


% %% Spectral Entropy
% LP=length(accP(:,1));
% fP=sr/LP*(0:(LP/2));
% 
% LF=length(accF(:,1));
% fF=sr/LF*(0:(LF/2));
% 
% % Spettro
% YP=fft(accP);
% xdftP=YP(1:LP/2+1,:);
% spettro_accP=(1/(sr*LP))*abs(xdftP).^2;
% spettro_accP(2:end-1,:)=2*spettro_accP(2:end-1,:);
% 
% YF=fft(accF);
% xdftF=YF(1:LF/2+1,:);
% spettro_accF=(1/(sr*LF))*abs(xdftF).^2;
% spettro_accF(2:end-1,:)=2*spettro_accF(2:end-1,:);
% 
% 
% pP=zeros(length(spettro_accP),1);
% spettro_accP_tot=sum(spettro_accP(:,1));
% 
% for i=1:length(spettro_accP)
%     pP(i)=spettro_accP(i,1)/spettro_accP_tot;
% end
% 
% entropiaP=-sum(pP.*log2(pP));
% 
% 
% pF=zeros(length(spettro_accF),1);
% spettro_accF_tot=sum(spettro_accF(:,1));
% 
% for i=1:length(spettro_accF)
%     pF(i)=spettro_accF(i,1)/spettro_accF_tot;
% end
% 
% entropiaF=-sum(pF.*log2(pF));
% 
% 
% % Spettro no Media
% YP_noMedia=fft(accP-accP_media);
% xdftP_noMedia=YP_noMedia(1:LP/2+1,:);
% spettro_accP_noMedia=(1/(sr*LP))*abs(xdftP_noMedia).^2;
% spettro_accP_noMedia(2:end-1,:)=2*spettro_accP_noMedia(2:end-1,:);
% 
% YF_noMedia=fft(accF-accF_media);
% xdftF_noMedia=YF_noMedia(1:LF/2+1,:);
% spettro_accF_noMedia=(1/(sr*LF))*abs(xdftF_noMedia).^2;
% spettro_accF_noMedia(2:end-1,:)=2*spettro_accF_noMedia(2:end-1,:);
% 
% 
% pP_noMedia=zeros(length(spettro_accP_noMedia),1);
% spettro_accP_noMedia_tot=sum(spettro_accP_noMedia(:,1));
% 
% for i=1:length(spettro_accP_noMedia)
%     pP_noMedia(i)=spettro_accP_noMedia(i,1)/spettro_accP_noMedia_tot;
% end
% 
% entropiaP_noMedia=-sum(pP_noMedia.*log2(pP_noMedia));
% 
% 
% pF_noMedia=zeros(length(spettro_accF_noMedia),1);
% spettro_accF_noMedia_tot=sum(spettro_accF_noMedia(:,1));
% 
% for i=1:length(spettro_accF_noMedia)
%     pF_noMedia(i)=spettro_accF_noMedia(i,1)/spettro_accF_noMedia_tot;
% end
% 
% entropiaF_noMedia=-sum(pF_noMedia.*log2(pF_noMedia));
% 
% 
% figure(Name="Spettro")
% subplot(2,1,1)
% plot(fP,spettro_accP(:,1),LineWidth=1,Color="b");
% title("Spettro Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,entropiaP*ones(length(spettro_accP(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Spectral Entropy",FontName=font)
% subplot(2,1,2)
% plot(fF,spettro_accF(:,1),LineWidth=1,Color="r")
% title("Trasformata Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,entropiaF*ones(length(spettro_accF(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Spectral Entropy",FontName=font)
% 
% 
% figure(Name="Spettro No Media")
% subplot(2,1,1)
% plot(fP,spettro_accP_noMedia(:,1),LineWidth=1,Color="b");
% title("Spettro No Media Acc X Piano",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fP,entropiaP_noMedia*ones(length(spettro_accP_noMedia(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Spectral Entropy",FontName=font)
% subplot(2,1,2)
% plot(fF,spettro_accF_noMedia(:,1),LineWidth=1,Color="r")
% title("Spettro No Media Acc X Forte",FontName=font)
% xlabel("Hz",FontName=font)
% ylabel("X''(f)",FontName=font)
% grid
% hold on
% plot(fF,entropiaF_noMedia*ones(length(spettro_accF_noMedia(:,1)),1),LineWidth=1,Color="black")
% legend("Trasformata","Spectral Entropy",FontName=font)


 

%%


% %% Velocità Angolare
% 
% limX=[floor(min(min(vangP(:,1)),min(vangF(:,1)))),ceil(max(max(vangP(:,1)),max(vangF(:,1))))];
% limY=[floor(min(min(vangP(:,2)),min(vangF(:,2)))),ceil(max(max(vangP(:,2)),max(vangF(:,2))))];
% limZ=[floor(min(min(vangP(:,3)),min(vangF(:,3)))),ceil(max(max(vangP(:,3)),max(vangF(:,3))))];
% 
% f=figure;
% subplot(3,2,1)
% plot(tP,vangP(:,1),LineWidth=1,Color="r")
% title("Velocità Angolare Piano",FontName=font)
% subtitle("Roll",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limX)
% grid
% subplot(3,2,3)
% plot(tP,vangP(:,2),LineWidth=1,Color="g")
% subtitle("Pitch",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limY)
% grid
% subplot(3,2,5)
% plot(tP,vangP(:,3),LineWidth=1,Color="b")
% subtitle("Yaw",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limZ)
% grid
% 
% subplot(3,2,2)
% plot(tF,vangF(:,1),LineWidth=1,Color="r")
% title("Velocità Angolare Forte",FontName=font)
% subtitle("Roll",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limX)
% grid
% subplot(3,2,4)
% plot(tF,vangF(:,2),LineWidth=1,Color="g")
% subtitle("Pitch",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limY)
% grid
% subplot(3,2,6)
% plot(tF,vangF(:,3),LineWidth=1,Color="b")
% subtitle("Yaw",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limZ)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_VangXYZ.png")
% end


% %% Varianza Velocità Angolare
% vangP_var=movvar(vangP,40);
% vangF_var=movvar(vangF,40);
% 
% 
% limX=[floor(min(min(vangP_var(:,1)),min(vangF_var(:,1)))),ceil(max(max(vangP_var(:,1)),max(vangF_var(:,1))))];
% limZ=[floor(min(min(vangP_var(:,3)),min(vangF_var(:,3)))),ceil(max(max(vangP_var(:,3)),max(vangF_var(:,3))))];
% 
% f=figure;
% subplot(2,2,1)
% plot(tP,vangP_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare Piano",FontName=font)
% subtitle("Roll",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,3)
% plot(tP,vangP_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limZ)
% grid
% 
% subplot(2,2,2)
% plot(tF,vangF_var(:,1),LineWidth=1,Color="r")
% title("Varianza Velocità Angolare Forte",FontName=font)
% subtitle("Roll",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limX)
% grid
% subplot(2,2,4)
% plot(tF,vangF_var(:,3),LineWidth=1,Color="b")
% subtitle("Yaw",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("deg/s",FontName=font)
% ylim(limZ)
% grid
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_varVangXZ.png")
% end


% %% Velocità
% velP=cumsum(accP(:,1));
% velF=cumsum(accF(:,1));
% 
% lim=[floor(min(min(velP),min(velF))),ceil(max(max(velP),max(velF)))];
% 
% f=figure;
% subplot(2,1,1)
% plot(tP,velP,LineWidth=1)
% title("Velocità",FontName=font)
% subtitle("Piano",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s",FontName=font)
% ylim(lim)
% grid
% hold on
% plot(tP,zeros(length(velP)),LineWidth=.1,Color="black")
% subplot(2,1,2)
% plot(tF,velF,LineWidth=1)
% subtitle("Forte",FontName=font)
% xlabel("t(s)",FontName=font)
% ylabel("m/s",FontName=font)
% ylim(lim)
% grid
% hold on
% plot(tF,zeros(length(velF)),LineWidth=.1,Color="black")
% 
% if(stampa)
%     exportgraphics(f,figurePath+"lungaFP_velX.png")
% end


%% Funzioni
function[sqm] = movrms(f,n)
newf=[zeros(n/2,3);f;zeros(n/2,3)];

sqm=zeros(length(f),3);

for i=1:length(f)
    sqm(i,:)=rms(newf(i:i+n,:));
end

end

function[kurt] = movkurt(f,n)
kurt=zeros(length(f),3);

for i=1:n
    kurt(i,1)=kurtosis([zeros(n-i,1);f(1:i,1)]);
    kurt(i,2)=kurtosis([zeros(n-i,1);f(1:i,2)]);
    kurt(i,3)=kurtosis([zeros(n-i,1);f(1:i,3)]);
end

for i=n+1:length(f)
    kurt(i,1)=kurtosis(f(i-n:i,1));
    kurt(i,2)=kurtosis(f(i-n:i,2));
    kurt(i,3)=kurtosis(f(i-n:i,3));
end

end

function[skew] = movskw(f,n)
skew=zeros(length(f),3);

for i=1:n
    skew(i,1)=skewness([zeros(n-i,1);f(1:i,1)]);
    skew(i,2)=skewness([zeros(n-i,1);f(1:i,2)]);
    skew(i,3)=skewness([zeros(n-i,1);f(1:i,3)]);
end

for i=n+1:length(f)
    skew(i,1)=skewness(f(i-n:i,1));
    skew(i,2)=skewness(f(i-n:i,2));
    skew(i,3)=skewness(f(i-n:i,3));
end

end


%% Grafici
function stampa(t,fun,tit,ax,xlbl,ylbl)
n=length(ax);
limY=ones(n,2);

textLimY=ones(n,5);

for i=1:n
    limY(i,:)=[floor(min(fun(:,i))), ceil(max(fun(:,i)))];
    textLimY(i,:)=(limY(i,2)-1)*ones(5,1);
end

stampa_gen(t,fun,tit,ax,xlbl,ylbl,limY,textLimY)
end

function stampa_gen(t,fun,tit,ax,xlbl,ylbl,limY,textLimY)
n=length(ax);
font="Times New Roman";

linee=ones(5,2);
linee(1,:)=[1 1];
linee(2,:)=[15 15];
linee(3,:)=[19.5 19.5];
linee(4,:)=[23.75 23.75];
linee(5,:)=[25.5 25.5];
% linee(6,:)=[29.75 29.75];

str={'Accelerate','Idle','Accelerate','Idle','Brake'};

for i=1:n
    if (ax(i)=='X')
        c(i)='r';
    elseif (ax(i)=='Y')
        c(i)='g';
    else
        c(i)='b';
    end
end

f=figure;
for i=1:n
    subplot(n,1,i)
    plot(t,fun(:,i),LineWidth=1,color=c(i))
    if(i==1)
        title(tit,FontName=font)
    end
    subtitle(ax(i),FontName=font)
    xlabel(xlbl,FontName=font)
    ylabel(ylbl,FontName=font)
    ylim(limY(i,:))
    grid
    hold on
    for j=1:length(linee)
        line(linee(j,:),limY(i,:), 'Color','black')
    end
    text(linee(:,1)+0.5*[1,1,1,1/2,1]',textLimY(i,:),str,FontSize=6.5, FontWeight="bold",FontName=font)

end

% exportgraphics(f,"D:\Users\Daniele\Desktop\"+tit+".png")

end
