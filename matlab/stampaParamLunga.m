clear all
close all
clc

path=".\dati\lunga_forte\";

rilievo=2;

[gzRot,gMedio] = GZRot(path);

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*gzRot*9.81/-gMedio;
vang_rotta=(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1);

[t,acc,vang,mag]=AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta);


%% Accelerazione
stampa(t,acc,"Accelerazione",['X','Y'],'t(s)','m/s^2')


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


% %% Deviazione Standard Accelerazione
% acc_std=movstd(acc,40);
%
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_std(:,i))),ceil(max(acc_std(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2);
%
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-0.5)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
%
% stampa_gen(t,acc_std(:,1:2),"AccStd",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Scarto Quadratico Medio Accelerazione
% acc_rms=movrms(acc,40);
%
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_rms(:,i))),ceil(max(acc_rms(:,i)))];
% end
% limY(1,2)=limY(1,2)+0.5;
% limY(2,2)=limY(2,2);
%
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
%
% stampa_gen(t,acc_rms(:,1:2),"AccRms",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Kurtosi Accelerazione
% acc_krt=movkurt(acc,40);
% stampa(t,acc_krt,"AccKurtosi",['X','Y'],'t(s)','m/s^2')


% %% Skewness Accelerazione
% acc_skw=movskw(acc,40);
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_skw(:,i))),ceil(max(acc_skw(:,i)))];
% end
% limY(1,2)=limY(1,2);
% limY(2,2)=limY(2,2)+0.5;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-0.5)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_skw(:,1:2),"AccSkewness",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Shape Factor Accelerazione
% acc_rms=movrms(acc,40);
% acc_arv=movmean(abs(acc),40);
% 
% acc_shf=acc_rms./acc_arv;
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_shf(:,i))),ceil(max(acc_shf(:,i)))];
% end
% limY(1,2)=limY(1,2)-0.5;
% limY(2,1)=limY(2,1)+0.5;
% limY(2,2)=limY(2,2)+0.5;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-0.5)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_shf(:,1:2),"AccShape",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Crest Factor Accelerazione
% acc_max=movmax(acc,40);
% acc_rms=movrms(acc,40);
% 
% acc_crf=acc_max./acc_rms;
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_crf(:,i))),ceil(max(acc_crf(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2)+1;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-0.5)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-0.5)*ones(5,1);
% 
% stampa_gen(t,acc_crf(:,1:2),"AccCrest",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Impulse Factor Accelerazione
% acc_max=movmax(acc,40);
% acc_arv=movmean(abs(acc),40);
% 
% acc_impf=acc_max./acc_arv;
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_impf(:,i))),ceil(max(acc_impf(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2)+1;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-1)*ones(5,1);
% 
% stampa_gen(t,acc_impf(:,1:2),"AccImp",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Margin Factor Accelerazione
% acc_max=movmax(acc,40);
% acc_arvq=movmean(sqrt(abs(acc)),40).^2;
% 
% acc_mrgf=acc_max./acc_arvq;
% 
% stampa(t,acc_mrgf(:,1:2),"AccMargin",['X','Y'],'t(s)','m/s^2')


% %% Max Accelerazione
% n_max=10;
% acc_max=movmax(acc,n_max);
% 
% stampa(t,acc_max,"AccMax",['X','Y'],'t(s)','m/s^2')


% %% Min Accelerazione
% n_min=10;
% acc_min=movmin(acc,n_min);
% 
% limY=ones(3,2);
% for i=1:3
%     limY(i,:)=[floor(min(acc_min(:,i))),ceil(max(acc_min(:,i)))];
% end
% limY(1,2)=limY(1,2)+1;
% limY(2,2)=limY(2,2)+1.5;
% 
% textLimY=ones(2,5);
% textLimY(1,:)=(limY(1,2)-1)*ones(5,1);
% textLimY(2,:)=(limY(2,2)-1)*ones(5,1);
% 
% stampa_gen(t,acc_min(:,1:2),"AccMin",['X','Y'],'t(s)','m/s^2',limY,textLimY)


% %% Peak Accelerazione
% n_Peak=10;
% acc_max=movmax(acc,n_Peak);
% acc_min=movmin(acc,n_Peak);
% 
% acc_peak=acc_max-acc_min;
% 
% stampa(t,acc_peak(:,1:2),"AccPeak",['X','Y'],'t(s)','m/s^2')


% %% Trasformata e Spettro Accelerazione
% stampa_freq(acc)


% %% Ampiezza Media Accelerazione
% stampa_AvAmp(acc)


% %% Frequency Centroid Accelerazione
% stampa_FrCen(acc)


% %% Frequency Variance Accelerazione
% stampa_FrVar(acc)


% %% Spectral Entropy Accelerazione
% stampa_SpEnt(acc)



%% QUI


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
newf=[zeros(n/2,3);f;zeros(n/2,3)];

kurt=zeros(length(f),3);

for i=1:length(f)
    kurt(i,:)=kurtosis(newf(i:i+n,:));
end

end

function[skew] = movskw(f,n)
newf=[zeros(n/2,3);f;zeros(n/2,3)];

skew=zeros(length(f),3);

for i=1:length(f)
    skew(i,:)=skewness(newf(i:i+n,:));
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
linee(3,:)=[19.4 19.4];
linee(4,:)=[23.76 23.76];
linee(5,:)=[25.4 25.4];
% linee(6,:)=[29.76 29.76];

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

function stampa_freq(fun)
font="Times New Roman";
sr = 25; %sample rate

L=length(fun);
freq=sr/L*(0:(L/2));

Y_noMedia=fft(fun-movmean(fun,40));
P2_noMedia=abs(Y_noMedia/L);
trasform=P2_noMedia(1:(L/2+1),:);
trasform(2:end-1,:)=2*trasform(2:end-1,:);

f=figure;
subplot(2,1,1)
plot(freq,trasform(:,1),LineWidth=1,Color="r")
title("Trasformata",FontName=font)
subtitle("X''",FontName=font)
xlabel('Hz',FontName=font)
ylabel("X''(Hz)",FontName=font)
grid
subplot(2,1,2)
plot(freq,trasform(:,2),LineWidth=1,Color="g")
subtitle("Y''",FontName=font)
xlabel('Hz',FontName=font)
ylabel("Y''(Hz)",FontName=font)
grid

% exportgraphics(f,"D:\Users\Daniele\Desktop\Trasformata.png")

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]

for i=1:length(sezione)
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    f=figure;
    subplot(2,1,1)
    plot(freq,trasform(:,1),LineWidth=1,Color="r")
    title("Trasformata "+str(i),FontName=font)
    subtitle("X''",FontName=font)
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    ylim([0 3])
    grid
    subplot(2,1,2)
    plot(freq,trasform(:,2),LineWidth=1,Color="g")
    subtitle("Y''",FontName=font)
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    ylim([0 1.5])
    grid

    % exportgraphics(f,"D:\Users\Daniele\Desktop\Tr"+str(i)+".png")
end

end

function stampa_AvAmp(fun)
font="Times New Roman";
sr = 25; %sample rate

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]';

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    trasform_mean=mean(trasform);

    plot(freq,trasform_mean(1)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    grid
    hold on
end
title("Average Amplitude X",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\AvgAmpX"+str(i)+".png")

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    trasform_mean=mean(trasform);

    plot(freq,trasform_mean(2)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    grid
    hold on
end
title("Average Amplitude Y",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\AvgAmpY"+str(i)+".png")
end

function stampa_FrCen(fun)
font="Times New Roman";
sr = 25; %sample rate

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]';

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    trasform_cent=mean(trasform.*freq');

    plot(freq,trasform_cent(1)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    grid
    hold on
end
title("Frequency Centroid X",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\FrCentX"+str(i)+".png")

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    trasform_cent=mean(trasform.*freq');

    plot(freq,trasform_cent(2)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    grid
    hold on
end
title("Frequency Centroid Y",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\FrCentY"+str(i)+".png")
end

function stampa_FrVar(fun)
font="Times New Roman";
sr = 25; %sample rate

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]';

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    centroid=mean(trasform.*freq');
    trasform_var=sum((freq'-centroid).*trasform(:,1))/sum(trasform(:,1));

    plot(freq,trasform_var(1)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    grid
    hold on
end
title("Frequency Variance X",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\FrVarX"+str(i)+".png")

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    P2_noMedia=abs(Y_noMedia/L);
    trasform=P2_noMedia(1:(L/2+1),:);
    trasform(2:end-1,:)=2*trasform(2:end-1,:);

    centroid=mean(trasform.*freq');
    trasform_var=sum((freq'-centroid).*trasform(:,1))/sum(trasform(:,1));

    plot(freq,trasform_var(2)*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    grid
    hold on
end
title("Frequency Variance Y",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\FrVarY"+str(i)+".png")
end

function stampa_SpEnt(fun)
font="Times New Roman";
sr = 25; %sample rate

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]';

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    xdftP=Y_noMedia(1:L/2+1,:);
    spettro=(1/(sr*L))*abs(xdftP).^2;
    spettro(2:end-1,:)=2*spettro(2:end-1,:);

    p=zeros(length(spettro),1);
    spettro_tot=sum(spettro(:,1));

    for j=1:length(spettro)
        p(j)=spettro(j,1)/spettro_tot;
    end
    
    entropia=-sum(p.*log2(p));

    plot(freq,entropia*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    grid
    hold on
end
title("Spectral Entropy X",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\SpEntX"+str(i)+".png")

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    xdftP=Y_noMedia(1:L/2+1,:);
    spettro=(1/(sr*L))*abs(xdftP).^2;
    spettro(2:end-1,:)=2*spettro(2:end-1,:);

    p=zeros(length(spettro),1);
    spettro_tot=sum(spettro(:,2));

    for j=1:length(spettro)
        p(j)=spettro(j,2)/spettro_tot;
    end
    
    entropia=-sum(p.*log2(p));

    plot(freq,entropia*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    grid
    hold on
end
title("Spectral Entropy Y",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\SpEntY"+str(i)+".png")
end



