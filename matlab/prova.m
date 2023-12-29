clear
close all
clc

%% Import dati

path="db\palazzago\";

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - frenata leggera
% 3 - frenata brusca con rimbalzo indietro
% 4 - inchiodata
% 5 - frenata brusca (3 frenate)
% 6 - prevalentemente discesa, al centro non ho pedalato
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
for i=2:length(t)
    intervalloT(i)=t(i)-t(i-1);
end
disp("tempo di campionamento minimo: "+num2str(min(intervalloT(2:end))));
disp("tempo di campionamento massimo: "+num2str(max(intervalloT(2:end))));

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

f = (0:length(acc)-1)*25/length(acc);


% acc=highpass(acc,0.01,25);
% Accelerazione
accf=fft(acc);

% % X
% % il picco che si vede tra i 0 e i 5Hz non sono sicuro sia effetto delle
% % pedalate perchè sulla strada dove ho preso queste misure c'era un dosso e
% % in altri rilievi, dove il dosso non c'era, non si vedeva. Seguiranno
% % ulteriori verifiche
% figure
% plot(f,abs(accf(:,1)),LineWidth=1,Color="r");
% title("trasformata discreta di fourier accelerazione in X");
% xlabel("f (Hz)");
% ylabel("|X''(f)|");
% 
% % Y
% % il picco visibile tra 0 e 2Hz potrebbe essere dovuto all'effetto delle
% % pedalate
% figure
% plot(f,abs(accf(:,2)),LineWidth=1,Color="g");
% title("trasformata discreta di fourier accelerazione in Y");
% xlabel("f (Hz)");
% ylabel("|Y''(f)|");
% 
% % Z
% figure
% plot(f,abs(accf(:,3)),LineWidth=1,Color="b");
% title("trasformata discreta di fourier accelerazione in Z");
% xlabel("f (Hz)");
% ylabel("|Z''(f)|");


% Velocità Angolare
vangf=fft(vang);

% % Roll (rotazione attorno a x)
% figure
% plot(f,abs(vangf(:,1)),LineWidth=1,Color="r");
% title("trasformata discreta di fourier Roll'");
% xlabel("f (Hz)");
% ylabel("|Roll'(f)|");
%
% % Pitch (rotazione attorno a y)
% figure
% plot(f,abs(vangf(:,2)),LineWidth=1,Color="g");
% title("trasformata discreta di fourier Pitch'");
% xlabel("f (Hz)");
% ylabel("|Pitch'(f)|");
%
% % Yaw (rotazione attorno a z)
% figure
% plot(f,abs(vangf(:,3)),LineWidth=1,Color="b");
% title("trasformata discreta di fourier Yaw'");
% xlabel("f (Hz)");
% ylabel("|Yaw'(f)|");

%% Filtraggio dati
sr=25; % frequenza di campionamento (sample rate) del sensore

lp1=1; % frequenza alla quale il filtro passa-basso esegue il taglio
hp1=0; % frequenza alla quale il filtro passa-alto esegue il taglio

filteredAcc1=lowpass(acc,lp1,sr);
% filteredAcc1=highpass(filteredAcc1,hp1,sr);
% multiPlotta3(t,acc,filteredAcc,"accelerazione","accelerazione filtrata");

filteredVang1=lowpass(vang,lp1,sr);
% filteredVang1=highpass(filteredVang1,hp1,sr);
% multiPlotta3(t,vang,filteredVang,"velocità angolare","velocità angolare filtrata");

plotta3(t,filteredAcc1,"accelerazione filtrata tra "+num2str(hp1)+" e "+num2str(lp1)+"Hz");
plotta3(t,filteredVang1,"velocità angolare filtrata tra "+num2str(hp1)+" e "+num2str(lp1)+"Hz");

moveAcc1=movmean(filteredAcc1,[15,0]);
plotta3(t,moveAcc1,"media mobile accelerazione");




lp2=5; % frequenza alla quale il filtro passa-basso esegue il taglio
hp2=1; % frequenza alla quale il filtro passa-alto esegue il taglio

filteredAcc2=lowpass(acc,lp2,sr);
filteredAcc2=highpass(filteredAcc2,hp2,sr);
% multiPlotta3(t,acc,filteredAcc,"accelerazione","accelerazione filtrata");

filteredVang2=lowpass(vang,lp2,sr);
filteredVang2=highpass(filteredVang2,hp2,sr);
% multiPlotta3(t,vang,filteredVang,"velocità angolare","velocità angolare filtrata");

plotta3(t,filteredAcc2,"accelerazione filtrata tra "+num2str(hp2)+" e "+num2str(lp2)+"Hz");
plotta3(t,filteredVang2,"velocità angolare filtrata tra "+num2str(hp2)+" e "+num2str(lp2)+"Hz");

moveAcc=movmean(filteredAcc2,[15,0]);
plotta3(t,moveAcc,"media mobile accelerazione");