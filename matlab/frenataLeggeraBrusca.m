clear
close all
clc

%% Import dati
path="db\palazzago\";
[gzRot,gMedio] = GZRot(path);

% 2 - frenata leggera
% 4 - inchiodata
% 5 - frenata brusca (3 frenate)
rLeggero=2; % rilievo frenata leggera (quindi più lunga)
rBrusco=4; % rilievo frenata brusca

dbl=importdata(path + "BlueCoin_Log_N00"+rLeggero+".csv").data;
dbb=importdata(path + "BlueCoin_Log_N00"+rBrusco+".csv").data;

inizioL=1;
fineL=length(dbl);

inizioB=1;
fineB=length(dbb);


tempoL=dbl(inizioL:fineL,1)*1e-3;
tempoL=tempoL-tempoL(1);

tempoB=dbb(inizioB:fineB,1)*1e-3;
tempoB=tempoB-tempoB(1);


% estrazione dati accelerometro (mg)
accelerazioneL=dbl(inizioL:fineL,2:4)*gzRot*9.81/-gMedio;
accelerazioneB=dbb(inizioB:fineB,2:4)*gzRot*9.81/-gMedio;

StampaAcc(tempoL,tempoB,accelerazioneL,accelerazioneB,"Accelerazione","Accelerazione Frenata Leggera","Accelerazione Frenata Brusca")


%% Accelerazione Filtrata
sr = 25;

% da 0 a 1Hz
lowFreq_lp=1;
lowFreq_hp=0;
lowFreqAccL=lowpass(accelerazioneL,lowFreq_lp,sr);
lowFreqAccB=lowpass(accelerazioneB,lowFreq_lp,sr);

if(lowFreq_hp ~= 0)
    lowFreqAccL=highpass(lowFreqAccL,lowFreq_hp,sr);
    lowFreqAccB=highpass(lowFreqAccB,lowFreq_hp,sr);
end

StampaAcc(tempoL,tempoB,lowFreqAccL,lowFreqAccB,"Accelerazione Frenata Leggera Filtrata tra "+num2str(lowFreq_hp)+" e "+num2str(lowFreq_lp)+"Hz","Accelerazione Frenata Leggera Filtrata","Accelerazione Frenata Brusca Filtrata")

n_val=25;
lowFreqAccL_media=movmean(lowFreqAccL,[n_val,0]);
lowFreqAccB_media=movmean(lowFreqAccB,[n_val,0]);

StampaAcc(tempoL,tempoB,lowFreqAccL_media,lowFreqAccB_media,"Media Mobile Accelerazione Frenata Leggera Filtrata tra "+num2str(lowFreq_hp)+" e "+num2str(lowFreq_lp)+"Hz","Accelerazione Frenata Leggera Filtrata","Accelerazione Frenata Brusca Filtrata")

% da 1 a 5Hz
highFreq_lp=5;
highFreq_hp=1;
highFreqAccL=lowpass(accelerazioneL,highFreq_lp,sr);
highFreqAccL=highpass(highFreqAccL,highFreq_hp,sr);
highFreqAccB=lowpass(accelerazioneB,highFreq_lp,sr);
highFreqAccB=highpass(highFreqAccB,highFreq_hp,sr);

StampaAcc(tempoL,tempoB,highFreqAccL,highFreqAccB,"Accelerazione Frenata Leggera Filtrata tra "+num2str(highFreq_hp)+" e "+num2str(highFreq_lp)+"Hz","Accelerazione Frenata Leggera Filtrata","Accelerazione Frenata Brusca Filtrata")


%% Velocità
velL=cumsum(accelerazioneL)*0.04;
velB=cumsum(accelerazioneB)*0.04;

StampaVel(tempoL,tempoB,velL,velB,"Velocità","Velocità Frenata Leggera","Velocità Frenata Brusca")


%% Velocità Angolare
vAngolareL=dbl(inizioL:fineL,5:7)*2*pi/360*1e-3;
vAngolareB=dbb(inizioB:fineB,5:7)*2*pi/360*1e-3;

StampaVang(tempoL,tempoB,vAngolareL,vAngolareB,"Velocità Angolare","Velocità Angolare Frenata Leggera","Velocità Angolare Frenata Brusca")


%% Integrale Velocità Angolare
angL=cumsum(vAngolareL)*0.04;
angB=cumsum(vAngolareB)*0.04;

StampaAng(tempoL,tempoB,angL,angB,"Integrale Velocità Angolare","Integrale Velocità Angolare Frenata Leggera","Integrale Velocità Angolare Frenata Brusca")


%% AHRS Filter
accelerazioneL=dbl(inizioL:fineL,2:4)*9.81/-gMedio;
accelerazioneB=dbb(inizioB:fineB,2:4)*9.81/-gMedio;

vAngolareL=deg2rad(dbl(inizioL:fineL,5:7)*1e-3);
vAngolareB=deg2rad(dbb(inizioB:fineB,5:7)*1e-3);

% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
magL=([dbl(inizioL:fineL,8),-dbl(inizioL:fineL,9),dbl(inizioL:fineL,10)]*1e-1);
magB=([dbb(inizioB:fineB,8),-dbb(inizioB:fineB,9),dbb(inizioB:fineB,10)]*1e-1);

fuse = ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');
% fuse.LinearAccelerationNoise=1e-1;
% fuse.GyroscopeDriftNoise=1e-1;
% reset(fuse);

[orientationL,angularVelocityL] = fuse(accelerazioneL,vAngolareL,magL);
[orientationB,angularVelocityB] = fuse(accelerazioneB,vAngolareB,magB);

StampaVang(tempoL,tempoB,angularVelocityL,angularVelocityB,"angularVelocity","angularVelocityP","angularVelocityF")
StampaAng(tempoL,tempoB,flip(rad2deg(unwrap(euler(orientationL,"ZYX","frame"))),2),flip(rad2deg(unwrap(euler(orientationB,"ZYX","frame"))),2),"Orientamento","Orientamento Leggero","Orientamento Brusco");


integrale_angularVelocityL=cumsum(angularVelocityL)*0.04;
integrale_angularVelocityB=cumsum(angularVelocityB)*0.04;

StampaAng(tempoL,tempoB,integrale_angularVelocityL,integrale_angularVelocityB,"Integrale angularVelocity","Integrale angularVelocity Leggera","Integrale angularVelocity Brusca")


%% Prisma che ruota
% orientationL
% figure
% pp=poseplot;
% timestamp = text(2,2,-2,num2str(tempoL(1)));
% xlabel('North')
% ylabel('East')
% for i=1:length(orientationL)
%     tic
%     % q=fuse(acc(i,:),vang(i,:),mag(i,:));
%     q = orientationL(i);
%     set(pp,"Orientation",q);
%     set(timestamp,"String",num2str(tempoL(i)))
%     drawnow
%     pause(1/(sr*2)-toc)
% end

% orientationB
% figure
% pp=poseplot;
% timestamp = text(2,2,-2,num2str(tempoB(1)));
% xlabel('North')
% ylabel('East')
% for i=1:length(orientationB)
%     tic
%     % q=fuse(acc(i,:),vang(i,:),mag(i,:));
%     q = orientationB(i);
%     set(pp,"Orientation",q);
%     set(timestamp,"String",num2str(tempoB(i)))
%     drawnow
%     pause(1/(sr*2)-toc)
% end




%% Funzioni
function Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,sottotitolo,etichettaX,etichettaY)

figure(Name=nome)
subplot(3,2,1)
plot(x1,y1(:,1),LineWidth=1,Color="r");
title(titolo1)
subtitle(sottotitolo(1))
grid
xlabel(etichettaX);
ylabel(etichettaY(1));

subplot(3,2,3)
plot(x1,y1(:,2),LineWidth=1,Color="g");
subtitle(sottotitolo(2))
grid
xlabel(etichettaX);
ylabel(etichettaY(2));

subplot(3,2,5)
plot(x1,y1(:,3),LineWidth=1,Color="b");
subtitle(sottotitolo(3))
grid
xlabel(etichettaX);
ylabel(etichettaY(3));


subplot(3,2,2)
plot(x2,y2(:,1),LineWidth=1,Color="r");
title(titolo2)
subtitle(sottotitolo(1))
grid
xlabel(etichettaX);
ylabel(etichettaY(1));

subplot(3,2,4)
plot(x2,y2(:,2),LineWidth=1,Color="g");
subtitle(sottotitolo(2))
grid
xlabel(etichettaX);
ylabel(etichettaY(2));

subplot(3,2,6)
plot(x2,y2(:,3),LineWidth=1,Color="b");
subtitle(sottotitolo(3))
grid
xlabel(etichettaX);
ylabel(etichettaY(3));

end

function StampaAcc(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"t(s)",["X''(m/s^2)","Y''(m/s^2)","Z''(m/s^2)"])
end

function StampaVel(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"t(s)",["X'(m/s)","Y'(m/s)","Z'(m/s)"])
end

function StampaPos(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"t(s)",["X(m)","Y(m)","Z(m)"])
end

function StampaVang(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"t(s)",["r'(rad/s)","p'(rad/s)","y'(rad/s)"])
end

function StampaAng(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"t(s)",["r(rad)","p(rad)","y(rad)"])
end

function StampaFreqAcc(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVAng(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end
