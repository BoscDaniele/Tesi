clear
close all
clc

%% Import dati
sr=25;
path="db\secchia\";
[gzRot,gMedio] = GZRot(path);

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - discesa via secchia
% 5 - curve
% 6 - pedalata tranquilla
rilievo=5;

% import dei dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db)-10*25;

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

%% Accelerazione
acc=db(inizio:fine,2:4)*gzRot*9.81/-gMedio;

low_lp=1;
low_hp=0;
lowFilteredAcc=lowpass(acc,low_lp,25);
if(low_hp ~= 0)
    lowFilteredAcc=highpass(lowFilteredAcc,low_hp,25);
end

lowFilteredAcc_mean=movmean(lowFilteredAcc,25);
StampaAcc(t,lowFilteredAcc_mean,"Accelerazione Media Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz","Accelerazione Media Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz")

high_lp=5;
high_hp=1;
highFilteredAcc=lowpass(acc,high_lp,25);
highFilteredAcc=highpass(highFilteredAcc,high_hp,25);
StampaAcc(t,highFilteredAcc,"Accelerazione Filtrata tra "+num2str(high_hp)+" e "+num2str(high_lp)+"Hz","Accelerazione Filtrata tra "+num2str(high_hp)+" e "+num2str(high_lp)+"Hz")


%% Velocità
vel=cumsum(acc)*0.04;
StampaVel(t,vel,"Velocità","Velocità")

%% Velocità Angolare
vang=deg2rad(db(inizio:fine,5:7)*1e-3);
StampaVang(t,vang,"Velocità Angolare","Velocità Angolare")

m_vang=movmean(vang,25);
StampaVang(t,m_vang,"Velocità Angolare Media","Velocità Angolare Media")


%% Angoli
ang=cumsum(vang)*0.04;
StampaPos(t,ang,"Integrale Velocità Angolare","Integrale Velocità Angolare")

m_ang=movmean(ang,25);
StampaAng(t,m_ang,"Integrale Velocità Angolare Media","Integrale Velocità Angolare Media")

% %% Prova DFT
% N=length(vang);
% vang_f=fft(vang);
% f=(0:N-1)*25/N;
% StampaFreqVAng(f,abs(vang_f.^2)/N,"dft vang","Trasformata Discreta di Fourier della Velocità Angolare")
% 
% StampaVang(t,vang,"vang","vang")
% 
% vang_filtrata=lowpass(vang,3,25);
% StampaVang(t,vang_filtrata,"vang_filtrata","vang_filtrata")
% 
% ang_filtrato=cumsum(vang_filtrata)*0.04;
% StampaAng(t,ang_filtrato,"ang_filtrato","ang_filtrato")



% %% Prova: elimino le variazioni troppo piccole della velocità angolare
% StampaAng(t,ang,"Angolo","Angolo")
% 
% diffAng=ang;
% newAng=ang;
% newAng(1,:)=[0,0,0];
% for i=2:length(vang)
%     diffAng(i,:)=(ang(i,:)+ang(i-1,:))/2;
% 
%     if(abs(diffAng(i,1)) <= 0.09)
%         newAng(i,1)=0;
%     end
% 
%     if(abs(diffAng(i,2)) <= 0.09)
%         newAng(i,2)=0;
%     end
% 
%     if(abs(diffAng(i,3)) <= 0.09)
%         newAng(i,3)=0;
%     end
% end
% StampaAng(t,newAng,"Angolo Nuovo","Angolo Nuovo")


%% AHRS Filter
acc=db(inizio:fine,2:4)*9.81/-gMedio;
vang=deg2rad(db(inizio:fine,5:7)*1e-3);

% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);

fuse = ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');
% fuse.LinearAccelerationNoise=1e-1;
% fuse.GyroscopeDriftNoise=1e-1;
% reset(fuse);

[orientation,angularVelocity] = fuse(acc,vang,mag);

StampaVang(t,angularVelocity,"angularVelocity","angularVelocity")
StampaAng(t,flip(rad2deg(unwrap(euler(orientation,"ZYX","frame"))),2),"Orientamento","Orientamento");


integrale_angularVelocity=cumsum(angularVelocity)*0.04;

StampaAng(t,integrale_angularVelocity,"Integrale angularVelocity","Integrale angularVelocity")


% %% Prisma che ruota
% figure
% pp=poseplot;
% timestamp = text(2,2,-2,num2str(t(1)));
% xlabel('North')
% ylabel('East')
% for i=1:length(orientation)
%     tic
%     % q=fuse(acc(i,:),vang(i,:),mag(i,:));
%     q = orientation(i);
%     set(pp,"Orientation",q);
%     set(timestamp,"String",num2str(t(i)))
%     drawnow
%     pause(1/(sr*2)-toc)
% end





%% Funzioni
function Stampa(x,y,nome,titolo,sottotitolo,etichettaX,etichettaY)
figure(Name=nome)
subplot(3,1,1)
plot(x,y(:,1),LineWidth=1,Color="r");
title(titolo)
subtitle(sottotitolo(1))
grid
xlabel(etichettaX);
ylabel(etichettaY(1));

subplot(3,1,2)
plot(x,y(:,2),LineWidth=1,Color="g");
subtitle(sottotitolo(2))
grid
xlabel(etichettaX);
ylabel(etichettaY(2));

subplot(3,1,3)
plot(x,y(:,3),LineWidth=1,Color="b");
subtitle(sottotitolo(3))
grid
xlabel(etichettaX);
ylabel(etichettaY(3));
end

function StampaAcc(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X''(m/s^2)","Y''(m/s^2)","Z''(m/s^2)"])
end

function StampaVel(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X'(m/s)","Y'(m/s)","Z'(m/s)"])
end

function StampaPos(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X(m)","Y(m)","Z(m)"])
end

function StampaVang(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"t(s)",["r'(rad/s)","p'(rad/s)","y'(rad/s)"])
end

function StampaAng(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"t(s)",["r(rad)","p(rad)","y(rad)"])
end

function StampaFreqAcc(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVang(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end
