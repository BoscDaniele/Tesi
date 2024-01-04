clear
close all
clc

%% Import dati
path="db\palazzago2\";

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - discesa + curve + alla fine salita
% 3 - frenata posteriore
% 4 - frenata anteriore
% 5 - percorso misto (leggera discesa + rotonda + leggera salita)
% 6 - discesa + molta salita
rilievo=5;

% import dei dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db);

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

% controllo il tempo di campionamento
% normalmente è di 0.04s ma è capitato che così non fosse
% for i=2:length(t)
%     intervalloT(i)=t(i)-t(i-1);
% end
% disp("tempo di campionamento minimo: "+num2str(min(intervalloT(2:end))));
% disp("tempo di campionamento massimo: "+num2str(max(intervalloT(2:end))));

% estrazione dati accelerometro (mg)
acc=db(inizio:fine,2:4);
% plotta3(t,acc,"accelerazioni");

% estrazione dati giroscopio e conversione in rad/s
vang=db(inizio:fine,5:7)*2*pi/360*1e-3;
% plotta3(t,vang,"velocità angolari");

%% Rotazione Sistema di Riferimento
% rotazione del sistema di riferimento del sensore al fine di farlo
% coincidere con quello della bicicletta

% calcolo matrice di rotazione attorno agli assi x, y, z

% gzRot contiene a matrice di rotazione

% gMedio contiene il valore della gravità moltiplicare il vettore
% accelerazione per 9.81/gMedio per convertire l'unità di misura da mg a
% m/s^2
[gzRot,gMedio] = GZRot(path);

% rotazione vettori
acc=acc*gzRot;
plotta3(t,acc,"accelerazioni ruotate");

vang=vang*gzRot; % non ne sono sicuro che funzioni così ma sicuramente in qualche modo vanno ruotati anche loro
plotta3(t,vang,"velocità angolari ruotate");

%% Trasformata di Fourier Discreta
% Per poter individuare quali sono le frequnze di interesse (e quindi quali
% filtrare) faccio la trasformata di fourier

f = (0:fine-1)*25/fine;

% Accelerazione
accf=highpass(acc,0.01,25);
accf=fft(accf);

% X
figure
plot(f,abs(accf(:,1)),LineWidth=1,Color="r");
title("trasformata discreta di fourier accelerazione in X");
xlabel("f (Hz)");
ylabel("|X''(f)|");

% Y
figure
plot(f,abs(accf(:,2)),LineWidth=1,Color="g");
title("trasformata discreta di fourier accelerazione in Y");
xlabel("f (Hz)");
ylabel("|Y''(f)|");

% Z
figure
plot(f,abs(accf(:,3)),LineWidth=1,Color="b");
title("trasformata discreta di fourier accelerazione in Z");
xlabel("f (Hz)");
ylabel("|Z''(f)|");


% Velocità Angolare
vangf=highpass(vang,0.01,25);
vangf=fft(vangf);

% Roll (rotazione attorno a x)
figure
plot(f,abs(vangf(:,1)),LineWidth=1,Color="r");
title("trasformata discreta di fourier Roll'");
xlabel("f (Hz)");
ylabel("|Roll'(f)|");

% Pitch (rotazione attorno a y)
figure
plot(f,abs(vangf(:,2)),LineWidth=1,Color="g");
title("trasformata discreta di fourier Pitch'");
xlabel("f (Hz)");
ylabel("|Pitch'(f)|");

% Yaw (rotazione attorno a z)
figure
plot(f,abs(vangf(:,3)),LineWidth=1,Color="b");
title("trasformata discreta di fourier Yaw'");
xlabel("f (Hz)");
ylabel("|Yaw'(f)|");

% multiPlotta3(f,abs(accf), abs(vangf), "trasformata accelerazione", "trasformata velocità angolare");

%% Filtraggio dati accelerazione
sr=25; % frequenza di campionamento (sample rate) del sensore
lpa=1; % frequenza alla quale il filtro passa-basso esegue il taglio
hpa=0; % frequenza alla quale il filtro passa-alto esegue il taglio

filteredAcc=lowpass(acc,lpa,sr);
% filteredAcc=highpass(filteredAcc,hpa,sr);

% multiPlotta3(t,acc,filteredAcc,"accelerazione","accelerazione filtrata");
plotta3(t,filteredAcc,"accelerazione filtrata tra "+num2str(hpa)+" e "+num2str(lpa)+"Hz");

moveAcc=movmean(filteredAcc,[25,0]);
plotta3(t,moveAcc,"media mobile accelerazione");



% filteredVang=lowpass(vang,lpa,sr);
% filteredVang=highpass(filteredVang,hpa,sr);
% multiPlotta3(t,vang,filteredVang,"velocità angolare","velocità angolare filtrata");
% plotta3(t,filteredVang,"velocità angolare filtrata tra "+num2str(hp)+" e "+num2str(lp)+"Hz");

