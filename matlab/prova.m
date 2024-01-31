clear
close all
clc

%% Import dati

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
fine=length(db);

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
StampaAcc(t,lowFilteredAcc,"Accelerazione Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz","Accelerazione Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz")
StampaAcc(t,movmean(acc,[40,0]),"acc mean","acc mean")


lowFilteredAcc_mean=movmean(lowFilteredAcc,25);
% StampaAcc(t,lowFilteredAcc_mean,"Accelerazione Media Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz","Accelerazione Media Filtrata tra "+num2str(low_hp)+" e "+num2str(low_lp)+"Hz")

high_lp=5;
high_hp=2;
highFilteredAcc=lowpass(acc,high_lp,25);
highFilteredAcc=highpass(highFilteredAcc,high_hp,25);
StampaAcc(t,highFilteredAcc,"Accelerazione Filtrata tra "+num2str(high_hp)+" e "+num2str(high_lp)+"Hz","Accelerazione Filtrata tra "+num2str(high_hp)+" e "+num2str(high_lp)+"Hz")


% %% Velocità
% vel=cumsum(acc)*0.04;
% StampaVel(t,vel,"Velocità","Velocità")


% %% Velocità Angolare
% vang=db(inizio:fine,5:7)*1e-3;
% StampaVang(t,vang,"Velocità Angolare","Velocità Angolare")
% 
% vang_f=lowpass(vang,0.5,25);
% StampaVang(t,vang_f,"Velocità Angolare Filtrata","Velocità Angolare Filtrata")
% 
% 
% m_vang=movmean(vang_f,25);
% % StampaVang(t,m_vang,"Velocità Angolare Media","Velocità Angolare Media")
% 
% 
% %% Angoli
% ang=cumsum(m_vang)*0.04;
% StampaAng(t,ang,"Posizione Angolare","Posizione Angolare")
% 
% m_ang=movmean(ang,25);
% StampaAng(t,m_ang,"Posizione Angolare Media","Posizione Angolare Media")


% %% Prova Integrale "progressivo"
% finestra=50;
% % Accelerazione
% newVel=zeros(length(acc),3);
% 
% for i = 1:finestra
%     newVel(finestra)=newVel(finestra)+acc(i);
% end
% 
% for i = finestra+1:length(acc)
%     newVel(i,1)=newVel(i-1,1)-acc(i-finestra,1)+acc(i,1);
%     newVel(i,2)=newVel(i-1,2)-acc(i-finestra,2)+acc(i,2);
%     newVel(i,3)=newVel(i-1,3)-acc(i-finestra,3)+acc(i,3);
% end
% 
% newVel=newVel.*0.04;
% StampaVel(t,newVel,"newVel","newVel")
% 
% 
% % Velocità angolare
% newAng=zeros(length(vang),3);
% 
% for i = 1:finestra
%     newAng(finestra)=newAng(finestra)+vang(i);
% end
% 
% for i = finestra+1:length(vang)
%     newAng(i,1)=newAng(i-1,1)-vang(i-finestra,1)+vang(i,1);
%     newAng(i,2)=newAng(i-1,2)-vang(i-finestra,2)+vang(i,2);
%     newAng(i,3)=newAng(i-1,3)-vang(i-finestra,3)+vang(i,3);
% end
% 
% newAng=newAng.*0.04;
% StampaAng(t,newAng,"newAng","newAng")


% %% Prova Soglia
% soglia=0.05;
% 
% dz=zeros(length(vang),3);
% 
% for i=2:length(vang_f)
%     dz(i,1)=(abs(vang(i,1))-abs(vang(i-1,1)));
%     dz(i,2)=(abs(vang(i,2))-abs(vang(i-1,2)));
%     dz(i,3)=(abs(vang(i,3))-abs(vang(i-1,3)));
% end
% 
% newVang=zeros(length(vang),3);
% 
% for i=1:length(dz)
%     if abs(dz(i,1))>soglia
%         newVang(i,1)=vang(i,1);
%     end
% 
%     if abs(dz(i,2))>soglia
%         newVang(i,2)=vang(i,2);
%     end
% 
%     if abs(dz(i,3))>soglia
%         newVang(i,3)=vang(i,3);
%     end
% end
% 
% StampaAng(t,newVang,"newVang","newVang");

%% Prova imufilter
acc=db(inizio:fine,2:4).*9.81/-gMedio;
% StampaAcc(t,acc,"acc","acc")

vang=deg2rad(db(inizio:fine,5:7).*1e-3);
% StampaVang(t,vang,"vang","vang")

% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);


% impostazioni
sr=25;

fuse = ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');
% fuse.LinearAccelerationNoise=1e-1;
% fuse.GyroscopeDriftNoise=1e-1;
% reset(fuse);

[orientation,angularVelocity] = fuse(acc,vang,mag);
StampaAng(t,flip(rad2deg(unwrap(euler(orientation,"ZYX","frame"))),2),"orientation","Orientation ahrsfilter")
% StampaVang(t,angularVelocity,"angularVelocity","Angular velocity with gyroscope bias removed in the sensor body coordinate system")

prova_orientation=euler(orientation,"XYZ","frame");

new_orientation=[prova_orientation(:,1:2),zeros(length(acc),1)];
no_zRot=quaternion(new_orientation,"euler","XYZ","frame");

% newAcc=zeros(length(acc),3);

% Sistema di riferimento IMU --> sistema di riferimento inerziale
newAcc=rotatepoint(orientation,acc);
newAcc = newAcc - [0,0,9.81];

newAng=euler(orientation,"ZYX","frame");
% rot=[-newAng(:,1),zeros(length(newAng),2)];
rot=-newAng;

newAcc3=rotatepoint(no_zRot,newAcc);
StampaAcc(t,newAcc3,"newAcc3","newAcc3")
g=[zeros(length(acc),2),ones(length(acc),1).*-9.81];

g_rot=rotatepoint(orientation,g);
% StampaAcc(t,g_rot,"g","g")


% newAcc3=acc+g_rot;

StampaAcc(t,acc,"acc","acc")
StampaAcc(t,newAcc,"newAcc ahrs","newAcc ahrs")
% StampaAcc(t,newAcc3,"newAcc ahrs ruotata","newAcc ahrs ruotata di -psi")


% % Trasformata Discreta di Fourier
% acc=highpass(acc,0.1,25);
% newAcc1=highpass(newAcc1,0.1,25);
% newAcc2=highpass(newAcc2,0.1,25);
% 
% acc_f=fft(acc);
% newAcc1_f=fft(newAcc1);
% newAcc2_f=fft(newAcc2);
% f=(0:length(t)-1)*25/length(t);
% 
% StampaFreqAcc(f,abs(acc_f),"frequenza acc","frequenza acc")
% StampaFreqAcc(f,abs(newAcc1_f),"frequenza imufilter","frequenza imufilter")
% StampaFreqAcc(f,abs(newAcc2_f),"frequenza ahrsfilter","frequenza ahrsfilter")


% % Velocità
% vel=cumsum(acc)*0.04;
newVel=cumsum(newAcc)*0.04;
newVel3=cumsum(newAcc3)*0.04;
% 
% StampaVel(t,vel,'vel','vel')
StampaVel(t,newVel,'newVel ahrs','newVel ahrs')
StampaVel(t,newVel3,'newVel ahrs routata','newVel ahrs ruotata')


% % Posizione
% pos=cumsum(vel);
% newPos1=cumsum(newVel1);
newPos=cumsum(newVel);
newPos3=cumsum(newVel3);
%
% StampaPos(t,pos,"pos","pos")
StampaPos(t,newPos,"newPos ahrs","newPos ahrs")
StampaPos(t,newPos3,"newPos ahrs ruotata","newPos ahrs ruotata")


%% Prisma che ruota

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
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"t(s)",["r(rad/s)","p(rad/s)","y(rad/s)"])
end

function StampaFreqAcc(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVang(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end
