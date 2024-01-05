clear
close all
clc

%% Import dati
path="db\secchia\";
[gzRot,gMedio] = GZRot(path);

rPiano=6; % rilievo pedalata piano
rForte=7; % rilievo pedalata forte

dbp=importdata(path + "BlueCoin_Log_N00"+rPiano+".csv").data;
dbf=importdata(path + "BlueCoin_Log_N00"+rForte+".csv").data;

% inizioP=20*25;
% fineP=length(dbp)-20*25;
%
% inizioF=10*25;
% fineF=length(dbf)-10*25;

inizioP=1;
fineP=length(dbp);

inizioF=1;
fineF=length(dbf);


tempoP=dbp(inizioP:fineP,1)*1e-3;
tempoP=tempoP-tempoP(1);

tempoF=dbf(inizioF:fineF,1)*1e-3;
tempoF=tempoF-tempoF(1);


% estrazione dati accelerometro (mg)
accelerazioneP=dbp(inizioP:fineP,2:4)*gzRot*9.81/-gMedio;
accelerazioneF=dbf(inizioF:fineF,2:4)*gzRot*9.81/-gMedio;

StampaAcc(tempoP,tempoF,accelerazioneP,accelerazioneF,"Accelerazione","Accelerazione Tranquila","Accelerazione Forte")


% estrazione dati giroscopio e conversione in rad/s
vAngolareP=dbp(inizioP:fineP,5:7)*gzRot*2*pi/360*1e-3;
vAngolareF=dbf(inizioF:fineF,5:7)*gzRot*2*pi/360*1e-3;

StampaVang(tempoP,tempoF,vAngolareP,vAngolareF,"Velocità Angolare","Velocità Angolare Tranquila","Velocità Angolare Forte")


%% Trasformata Discreta di Fourier Accelerazione
freqPiano=(0:length(tempoP)-1)*25/length(tempoP);
freqForte=(0:length(tempoF)-1)*25/length(tempoF);

TrasfAccPiano=fft(accelerazioneP);
TrasfAccForte=fft(accelerazioneF);

StampaFreqAcc(freqPiano,freqForte,abs(TrasfAccPiano),abs(TrasfAccForte),"Trasformata Accelerazione","Trasformata Accelerazione Tranquila","Trasformata Accelerazione Forte")


%% Accelerazione Filtrata
sr = 25;

% da 0 a 1Hz
lowFreq_lp=1;
lowFreq_hp=0;
lowFreqAccP=lowpass(accelerazioneP,lowFreq_lp,sr);
lowFreqAccF=lowpass(accelerazioneF,lowFreq_lp,sr);

if(lowFreq_hp ~= 0)
    lowFreqAccP=highpass(lowFreqAccP,lowFreq_hp,sr);
    lowFreqAccF=highpass(lowFreqAccF,lowFreq_hp,sr);
end

StampaAcc(tempoP,tempoF,lowFreqAccP,lowFreqAccF,"Accelerazione Filtrata tra "+num2str(lowFreq_hp)+" e "+num2str(lowFreq_lp)+"Hz","Accelerazione Tranquila Filtrata","Accelerazione Forte Filtrata")

n_val=25;
lowFreqAccP_media=movmean(lowFreqAccP,[n_val,0]);
lowFreqAccF_media=movmean(lowFreqAccF,[n_val,0]);

StampaAcc(tempoP,tempoF,lowFreqAccP_media,lowFreqAccF_media,"Media Mobile Accelerazione Filtrata tra "+num2str(lowFreq_hp)+" e "+num2str(lowFreq_lp)+"Hz","Accelerazione Tranquila Filtrata","Accelerazione Forte Filtrata")

% da 1 a 5Hz
highFreq_lp=5;
highFreq_hp=1;
highFreqAccP=lowpass(accelerazioneP,highFreq_lp,sr);
highFreqAccP=highpass(highFreqAccP,highFreq_hp,sr);
highFreqAccF=lowpass(accelerazioneF,highFreq_lp,sr);
highFreqAccF=highpass(highFreqAccF,highFreq_hp,sr);

StampaAcc(tempoP,tempoF,highFreqAccP,highFreqAccF,"Accelerazione Filtrata tra "+num2str(highFreq_hp)+" e "+num2str(highFreq_lp)+"Hz","Accelerazione Tranquila Filtrata","Accelerazione Forte Filtrata")


%% Velocità
velP=cumsum(accelerazioneP)*0.04;
velF=cumsum(accelerazioneF)*0.04;

StampaVel(tempoP,tempoF,velP,velF,"Velocità","Velocità Tranquila","Velocità Forte")

%% Posizione
posP=cumsum(velP)*0.04;
posF=cumsum(velF)*0.04;

StampaPos(tempoP,tempoF,posP,posF,"Posizione","Posizione Tranquila","Posizione Forte")


%% Trasformata Discreta di Fourier Velocità Angolare
TrasfVAngPiano=fft(vAngolareP);
TrasfVAngForte=fft(vAngolareF);

StampaFreqVAng(freqPiano,freqForte,abs(TrasfVAngPiano),abs(TrasfVAngForte),"Trasformata Velocità Angolare","Trasformata Velocità Angolare Tranquila","Trasformata Velocità Angolare Forte")


%% Velocità Angolare Filtrata
sr = 25;

% da 0 a 1Hz
lowFreq_lp=0.5;
lowFreq_hp=0;
lowFreqVAngP=lowpass(vAngolareP,lowFreq_lp,sr);
lowFreqVAngF=lowpass(vAngolareF,lowFreq_lp,sr);

if(lowFreq_hp ~= 0)
    lowFreqVAngP=highpass(lowFreqVAngP,lowFreq_hp,sr);
    lowFreqVAngF=highpass(lowFreqVAngF,lowFreq_hp,sr);
end

% StampaVang(tempoP,tempoF,lowFreqVAngP,lowFreqVAngF,"Velocità Angolare Filtrata tra "+num2str(lowFreq_hp)+" e "+num2str(lowFreq_lp)+"Hz","Velocità Angolare Tranquilla","Velocità Angolare Forte")


% da 1 a 5Hz
highFreq_lp=5;
highFreq_hp=1;
highFreqVAngP=lowpass(vAngolareP,highFreq_lp,sr);
highFreqVAngF=lowpass(vAngolareF,highFreq_lp,sr);

% StampaVang(tempoP,tempoF,highFreqVAngP,highFreqVAngF,"Velocità Angolare Filtrata tra "+num2str(highFreq_hp)+" e "+num2str(highFreq_lp)+"Hz","Velocità Angolare Tranquilla","Velocità Angolare Forte")

%% Posizione Angolare
angP=cumsum(vAngolareP)*0.04;
angF=cumsum(vAngolareF)*0.04;

StampaAng(tempoP,tempoF,angP,angF,"Posizione Angolare","Posizione Angolare Tranquilla","Posizione Angolare Forte")





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
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"t(s)",["r(rad/s)","p(rad/s)","y(rad/s)"])
end

function StampaFreqAcc(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVAng(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end
