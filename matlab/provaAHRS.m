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

sr = 25; %sample rate

%% Import Dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db);

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*9.81/-gMedio;
StampaAcc(t,acc*gzRot,"acc","acc")

vang=deg2rad(db(inizio:fine,5:7)*1e-3);
StampaVang(t,vang,"vang","vang")
ang=cumsum(vang)*0.04;
StampaAng(t,ang,"ang","ang")


% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);


%% ahrsFilter
fuse = ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');
% fuse.LinearAccelerationNoise=1e-1;
% fuse.GyroscopeDriftNoise=1e-1;
% reset(fuse);

[orientation,angularVelocity] = fuse(acc,vang,mag);
StampaAng(t,flip(rad2deg(unwrap(euler(orientation,"ZYX","frame"))),2),"orientation","Orientation ahrsfilter frame")
% StampaVang(t,angularVelocity,"angularVelocity","Angular velocity with gyroscope bias removed in the sensor body coordinate system")
% StampaAng(t,rad2deg(unwrap(euler(orientation,"XYZ","frame"))),"orientation","Orientation ahrsfilter prova XYZ")


StampaVang(t,angularVelocity,"angularVelocity","angularVelocity")
angular=cumsum(angularVelocity)*0.04;
StampaAng(t,angular,"angular","angular")


% Sistema di riferimento IMU --> sistema di riferimento inerziale
% newAcc=rotatepoint(orientation,acc);
% StampaAcc(t,newAcc,"newAcc ahrs","newAcc ahrs")
% 
% newAccFrame=rotateframe(orientation,acc);
% StampaAcc(t,newAccFrame,"newAccFrame ahrs","newAccFrame ahrs")


% provacon matrice di rotazione
fuse2 = ahrsfilter('SampleRate',sr,'OrientationFormat','Rotation matrix', ...
    'ReferenceFrame','NED');

[orientationMatrix,~] = fuse2(acc,vang,mag);


newG=zeros(length(acc),3);
for i=1:length(newG)
    newG(i,:)=[0,0,9.81]*orientationMatrix(:,:,i)';
end

StampaAcc(t,newG,"newG","newG")
accNoG=acc-newG;
StampaAcc(t,accNoG,"acc-g","acc-g")

accMedio=movmean(acc,[40,0]);
accNoGMedio=movmean(accNoG,[40,0]);
StampaAcc(t,accMedio*gzRot,"Accelerazione Media","Accelerazione Media")
StampaAcc(t,accNoGMedio*gzRot,"Accelerazione Media senza g","Accelerazione Media senza g")


% %% complementaryFilter
% 
% compfuse = complementaryFilter('SampleRate',sr,'OrientationFormat','quaternion','ReferenceFrame','NED');
% [compOrientation,compAngularVelocity] = compfuse(acc,vang,mag);
% 
% StampaAng(t,flip(rad2deg(unwrap(euler(compOrientation,"ZYX","frame"))),2),"compOrientation","Orientation Complementary Filter")
% StampaAng(t,rad2deg(unwrap(euler(compOrientation,"XYZ","frame"))),"compOrientation","Orientation Complementary Filter XYZ")



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
