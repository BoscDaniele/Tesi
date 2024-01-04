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

inizioP=20*25;
fineP=length(dbp)-20*25;

inizioF=10*25;
fineF=length(dbf)-10*25;

% inizioP=1;
% fineP=length(dbp);
% 
% inizioF=1;
% fineF=length(dbf);


tempoP=dbp(inizioP:fineP)*1e-3;
tempoP=tempoP-tempoP(1);

tempoF=dbf(inizioF:fineF)*1e-3;
tempoF=tempoF-tempoF(1);


% estrazione dati accelerometro (mg)
accelerazioneP=dbp(inizioP:fineP,2:4)*gzRot;
accelerazioneF=dbf(inizioF:fineF,2:4)*gzRot;

StampaAcc(tempoP,tempoF,accelerazioneP,accelerazioneF,"Accelerazione","Accelerazione Tranquila","Accelerazione Forte")


% estrazione dati giroscopio e conversione in rad/s
vAngolareP=dbp(inizioP:fineP,5:7)*2*pi/360*1e-3;
vAngolareF=dbf(inizioF:fineF,5:7)*2*pi/360*1e-3;

% StampaVang(tempoP,tempoF,vAngolareP,vAngolareF,"Velocità Angolare","Velocità Angolare Tranquila","Velocità Angolare Forte")


%% Trasformata Discreta di Fourier
% Trasformata Accelerazione
freqPiano=(0:length(tempoP)-1)*25/length(tempoP);
freqForte=(0:length(tempoF)-1)*25/length(tempoF);

TrasfAccPiano=fft(accelerazioneP);
TrasfAccForte=fft(accelerazioneF);

StampaFreqAcc(freqPiano,freqForte,abs(TrasfAccPiano),abs(TrasfAccForte),"Trasformata Accelerazione","Trasformata Accelerazione Tranquila","Trasformata Accelerazione Forte")


% Trasformmata Velocità Angolare
TrasfVAngPiano=fft(vAngolareP);
TrasfVAngForte=fft(vAngolareF);

% StampaFreqVAng(freqPiano,freqForte,abs(TrasfVAngPiano),abs(TrasfVAngForte),"Trasformata Velocità Angolare","Trasformata Velocità Angolare Tranquila","Trasformata Velocità Angolare Forte")


%% Accelerazione Filtrata
sr = 25;

% da 0 a 1Hz
lowFreq_lp=1;
lowFreqAccP=lowpass(accelerazioneP,lowFreq_lp,sr);
lowFreqAccF=lowpass(accelerazioneF,lowFreq_lp,sr);

StampaAcc(tempoP,tempoF,lowFreqAccP,lowFreqAccF,"Accelerazione Filtrata","Accelerazione Tranquila Filtrata","Accelerazione Forte Filtrata")



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
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"t(s)",["acc(mg)","acc(mg)","acc(mg)"])
end

function StampaVang(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"t(s)",["r'(rad/s)","p'(rad/s)","y'(rad/s)"])
end

function StampaFreqAcc(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVAng(x1,x2,y1,y2,nome,titolo1,titolo2)
Stampa(x1,x2,y1,y2,nome,titolo1,titolo2,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end
