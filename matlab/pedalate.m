clear
close all
clc

%% Import dati

% selezionare il database dalla quale caricare i dati
path="db\pedalate\";

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
rilievo=6;

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
lp=5; % frequenza alla quale il filtro passa-basso esegue il taglio
hp=1; % frequenza alla quale il filtro passa-alto esegue il taglio

filteredAcc=lowpass(acc,lp,sr);
% filteredAcc=highpass(filteredAcc,hp,sr);
% multiPlotta3(t,acc,filteredAcc,"accelerazione","accelerazione filtrata");

filteredVang=lowpass(vang,lp,sr);
filteredVang=highpass(filteredVang,hp,sr);
% multiPlotta3(t,vang,filteredVang,"velocità angolare","velocità angolare filtrata");

plotta3(t,filteredAcc,"accelerazione filtrata tra "+num2str(hp)+" e "+num2str(lp)+"Hz");
% plotta3(t,filteredVang,"velocità angolare filtrata tra "+num2str(hp)+" e "+num2str(lp)+"Hz");


% Nel caso l'ampiezza e la frequnza delle oscillazioni dovesse rivelarsi
% utile (come penso) nel determinare la cadenza e "l'impegno" del ciclista
% nel pedalare.
% La cadenza potrebbe essere determinata valutando quando
% l'accelerazione (una volta rimossa la componente non oscillante) passa
% dall'essere positiva all'essere negativa.
% L'impegno del ciclista (ampiezza delle oscillazioni) potrebbe essere
% valutato tramite l'inviluppo (sequenza dei massimi)? Sarà possibile, nel
% caso, farlo anche in real time?


%% Cadenza
cadenceX=zeros(length(acc),1);
cadenceY=zeros(length(acc),1);
tx1=1;
tx2=0;
ty1=1;
ty2=0;

for i=2:fine
    if(filteredAcc(i-1,1)>=0 && filteredAcc(i,1)<0)
        tx2=i;
        cadenceX(i)=1/((tx2-tx1)*0.04);
        tx1=i;
        % disp("rotazione pedale completata in "+num2str(1/cadenceX(i))+"s");
    end

    if(filteredAcc(i-1,2)>=0 && filteredAcc(i,2)<0)
        ty2=i;
        cadenceY(i)=1/((ty2-ty1)*0.04);
        ty1=i;
        % disp("rotazione pedale completata in "+num2str(1/cadenceY(i))+"s");
    end
end

for i=fine-1:-1:1
    if(cadenceX(i)==0)
        cadenceX(i)=cadenceX(i+1);
    end

    % if(cadenceY(i)==0)
    %     cadenceY(i)=cadenceY(i+1);
    % end
end

figure
plot(t,cadenceX,LineWidth=1,Color="g");
title("cadenza")
hold on
plot(t,filteredAcc(:,1)*9.81/-gMedio,LineWidth=1,Color="r");
grid
xlabel("t(s)");
ylabel("rotazioni/s - m/s^2");
legend("cadenza (rotazioni/s)", "accelerazione (m/s^2)")
