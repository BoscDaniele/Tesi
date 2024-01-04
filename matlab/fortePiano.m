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


tempoP=dbp(inizioP:fineP)*1e-3;
tempoP=tempoP-tempoP(1);

tempoF=dbf(inizioF:fineF)*1e-3;
tempoF=tempoF-tempoF(1);


% estrazione dati accelerometro (mg)
accelerazioneP=dbp(inizioP:fineP,2:4)*gzRot;
accelerazioneF=dbf(inizioF:fineF,2:4)*gzRot;

figure
subplot(3,2,1)
plot(tempoP,accelerazioneP(:,1),LineWidth=1,Color="r");
title("Accelerazione Tranquila")
subtitle("X")
grid
xlabel("t(s)");
ylabel("acc(mg)");

subplot(3,2,3)
plot(tempoP,accelerazioneP(:,2),LineWidth=1,Color="g");
subtitle("Y")
grid
xlabel("t(s)");
ylabel("acc(mg)");

subplot(3,2,5)
plot(tempoP,accelerazioneP(:,3),LineWidth=1,Color="b");
subtitle("Z")
grid
xlabel("t(s)");
ylabel("acc(mg)");


subplot(3,2,2)
plot(tempoF,accelerazioneF(:,1),LineWidth=1,Color="r");
title("Accelerazione Forte")
subtitle("X")
grid
xlabel("t(s)");
ylabel("acc(mg)");

subplot(3,2,4)
plot(tempoF,accelerazioneF(:,2),LineWidth=1,Color="g");
subtitle("Y")
grid
xlabel("t(s)");
ylabel("acc(mg)");

subplot(3,2,6)
plot(tempoF,accelerazioneF(:,3),LineWidth=1,Color="b");
subtitle("Z")
grid
xlabel("t(s)");
ylabel("acc(mg)");

% estrazione dati giroscopio e conversione in rad/s
vAngolareP=dbp(inizioP:fineP,5:7)*2*pi/360*1e-3;
vAngolareF=dbf(inizioF:fineF,5:7)*2*pi/360*1e-3;

% figure
% subplot(3,2,1)
% plot(tempoP,vAngolareP(:,1),LineWidth=1,Color="r");
% title("Velocità Angolare Tranquila")
% subtitle("Roll")
% grid
% xlabel("t(s)");
% ylabel("r'(rad/s)");
% 
% subplot(3,2,3)
% plot(tempoP,vAngolareP(:,2),LineWidth=1,Color="g");
% subtitle("Pitch")
% grid
% xlabel("t(s)");
% ylabel("p'(rad/s)");
% 
% subplot(3,2,5)
% plot(tempoP,vAngolareP(:,3),LineWidth=1,Color="b");
% subtitle("Yaw")
% grid
% xlabel("t(s)");
% ylabel("y'(rad/s)");
% 
% 
% subplot(3,2,2)
% plot(tempoF,vAngolareF(:,1),LineWidth=1,Color="r");
% title("Velocità Angolare Forte")
% subtitle("Roll")
% grid
% xlabel("t(s)");
% ylabel("r'(rad/s)");
% 
% subplot(3,2,4)
% plot(tempoF,vAngolareF(:,2),LineWidth=1,Color="g");
% subtitle("Pitch")
% grid
% xlabel("t(s)");
% ylabel("p'(rad/s)");
% 
% subplot(3,2,6)
% plot(tempoF,vAngolareF(:,3),LineWidth=1,Color="b");
% subtitle("Yaw")
% grid
% xlabel("t(s)");
% ylabel("y'(rad/s)");

%% Trasformata Discreta di Fourier
% Trasformata Accelerazione
freqPiano=(0:length(tempoP)-1)*25/length(tempoP);
freqForte=(0:length(tempoF)-1)*25/length(tempoF);

TrasfAccPiano=fft(accelerazioneP);
TrasfAccForte=fft(accelerazioneF);

figure
subplot(3,2,1)
plot(freqPiano,abs(TrasfAccPiano(:,1)),LineWidth=1,Color="r");
title("Trasformata Accelerazione Tranquila")
subtitle("X")
grid
xlabel("f(Hz)");
ylabel("|X''(f)|");

subplot(3,2,3)
plot(freqPiano,abs(TrasfAccPiano(:,2)),LineWidth=1,Color="g");
subtitle("Y")
grid
xlabel("f(Hz)");
ylabel("|Y''(f)|");

subplot(3,2,5)
plot(freqPiano,abs(TrasfAccPiano(:,3)),LineWidth=1,Color="b");
subtitle("Z")
grid
xlabel("f(Hz)");
ylabel("|Z''(f)|");


subplot(3,2,2)
plot(freqForte,abs(TrasfAccForte(:,1)),LineWidth=1,Color="r");
title("Trasformata Accelerazione Forte")
subtitle("X")
grid
xlabel("f(Hz)");
ylabel("|X''(f)|");

subplot(3,2,4)
plot(freqForte,abs(TrasfAccForte(:,2)),LineWidth=1,Color="g");
subtitle("Y")
grid
xlabel("f(Hz)");
ylabel("|Y''(f)|");

subplot(3,2,6)
plot(freqForte,abs(TrasfAccForte(:,3)),LineWidth=1,Color="b");
subtitle("Z")
grid
xlabel("f(Hz)");
ylabel("|Z''(f)|");


% Trasformmata Velocità Angolare
TrasfVAngPiano=fft(vAngolareP);
TrasfVAngForte=fft(vAngolareF);

figure
subplot(3,2,1)
plot(freqPiano,abs(TrasfVAngPiano(:,1)),LineWidth=1,Color="r");
title("Trasformata Accelerazione Tranquila")
subtitle("Roll")
grid
xlabel("f(Hz)");
ylabel("|r'(f)|");

subplot(3,2,3)
plot(freqPiano,abs(TrasfVAngPiano(:,2)),LineWidth=1,Color="g");
subtitle("Pitch")
grid
xlabel("f(Hz)");
ylabel("|p'(f)|");

subplot(3,2,5)
plot(freqPiano,abs(TrasfVAngPiano(:,3)),LineWidth=1,Color="b");
subtitle("Yaw")
grid
xlabel("f(Hz)");
ylabel("|y'(f)|");


subplot(3,2,2)
plot(freqForte,abs(TrasfVAngForte(:,1)),LineWidth=1,Color="r");
title("Trasformata Accelerazione Forte")
subtitle("Roll")
grid
xlabel("f(Hz)");
ylabel("|r'(f)|");

subplot(3,2,4)
plot(freqForte,abs(TrasfVAngForte(:,2)),LineWidth=1,Color="g");
subtitle("Pitch")
grid
xlabel("f(Hz)");
ylabel("|p'(f)|");

subplot(3,2,6)
plot(freqForte,abs(TrasfVAngForte(:,3)),LineWidth=1,Color="b");
subtitle("Yaw")
grid
xlabel("f(Hz)");
ylabel("|y'(f)|");





