clear
close all
clc

%% Import dati

% selezionare il database dalla quale caricare i dati
path="dbdm\pedalate\";

% in questa cartella ci sono 7 file (numerati da 0 a 6), i primi due
% servono solo a orientare il sensore, negli altri sfortunatamente non
% ricordo cosa ho fatto nello specifico (non ricordo in che ordine), ma
% queste misure erano state prese con lo scolpo di verificare se in base
% alla velocità cambiava la frequenza e l'ampiezza delle accelerazioni.
% Sono state eseguite anche delle prove per vedere se il cambio dei
% rapporti ("marcia") risultava visibile all'interno del grafico delle
% accelerazioni (variazione frequenza e picco nell'ampiezza delle
% oscillazioni). Purtroppo la presenza di un dosso rende le misure meno
% affidabili.

% selezionare il rilievo da caricare
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
% plotta3(t,acc,"accelerazioni ruotate");

vang=vang*gzRot; % non ne sono sicuro che funzioni così ma sicuramente in qualche modo vanno ruotati anche loro
% plotta3(t,vang,"velocità angolari ruotate");

%% Trasformata di Fourier Discreta
% Per poter individuare quali sono le frequnze di interesse (e quindi quali
% filtrare) faccio la trasformata di fourier

f = (0:length(acc)-1)*25/length(acc);

% Accelerazione
accf=fft(acc);

% X
% il picco che si vede tra i 0 e i 5Hz non sono sicuro sia effetto delle
% pedalate perchè sulla strada dove ho preso queste misure c'era un dosso e
% in altri rilievi, dove il dosso non c'era, non si vedeva. Seguiranno
% ulteriori verifiche
figure
plot(f,abs(accf(:,1)),LineWidth=1,Color="r");
title("trasformata discreta di fourier accelerazione in X");
xlabel("f (Hz)");
ylabel("|X''(f)|");

% Y
% il picco visibile tra 0 e 2Hz potrebbe essere dovuto all'effetto delle
% pedalate
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
vangf=fft(vang);

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

%% Filtraggio dati
% lp=5; % frequenza alla quale il filtro passabasso esegue il taglio
% 
% acc=lowpass(acc,lp,25);
% vang=lowpass(vang,lp,25);
% 
% plotta3(t,acc,"accelerazione filtrata a "+num2str(lp)+"Hz");
% plotta3(t,vang,"velocità angolare filtrata a "+num2str(lp)+"Hz");


% Nel caso l'ampiezza e la frequnza delle oscillazioni dovesse rivelarsi
% utile (come penso) nel determinare la cadenza e "l'impegno" del ciclista
% nel pedalare la cadenza potrebbe essere determinata valutando quando
% l'accelerazione (una volta rimossa la componente non oscillante) passa
% dall'essere positiva all'essere negativa.
% L'impegno del ciclista (ampiezza delle oscillazioni) potrebbe essere
% valutato tramite l'inviluppo (sequenza dei massimi)? Sarà possibile, nel
% caso, farlo anche in real time?


