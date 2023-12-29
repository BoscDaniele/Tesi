clear
close all
clc

%% Import dati

path="dbdm\palazzago\";

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - frenata leggera
% 3 - frenata brusca con rimbalzo indietro
% 4 - inchiodata
% 5 - frenata brusca (3 frenate)
% 6 - prevalentemente discesa, al centro non ho pedalato
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
lp=5; % frequenza alla quale il filtro passa-basso esegue il taglio
hp=0.01; % frequenza alla quale il filtro passa-alto esegue il taglio

filteredAcc=lowpass(acc,lp,sr);
filteredAcc=highpass(filteredAcc,hp,sr);
% multiPlotta3(t,acc,filteredAcc,"accelerazione","accelerazione filtrata");

filteredVang=lowpass(vang,lp,sr);
filteredVang=highpass(filteredVang,hp,sr);
% multiPlotta3(t,vang,filteredVang,"velocità angolare","velocità angolare filtrata");

plotta3(t,filteredAcc,"accelerazione filtrata tra "+num2str(hp)+" e "+num2str(lp)+"Hz");
plotta3(t,filteredVang,"velocità angolare filtrata tra "+num2str(hp)+" e "+num2str(lp)+"Hz");

moveAcc=movmean(filteredAcc,[15,0]);
plotta3(t,moveAcc,"media mobile accelerazione");




% intervallo=30;
% 
% for i=1:ceil(fine/(25*intervallo))
%     start=(i-1)*(25*intervallo)+1;
%     if(i*(25*intervallo)-1>fine)
%         final=fine;
%     else
%         final=i*(25*intervallo);
%     end
% 
%    multiPlotta3(t(start:final),filteredAcc(start:final,:),filteredVang(start:final,:),"accelerazione da "+num2str(t(start))+"s a "+num2str(t(final))+"s","velocità angolare da "+num2str(t(start))+"s a "+num2str(t(final))+"s")
% 
% end



vel=cumsum(acc)*0.04;
pos=cumsum(vel)*0.04;

filteredVel=highpass(vel,hp,sr);
filteredVel=lowpass(filteredVel,lp,sr);

filteredPos=highpass(pos,hp,sr);
filteredPos=lowpass(filteredPos,lp,sr);

% plotta3(t,filteredVel,"filteredVel");
% plotta3(t,filteredPos,"filteredPos");


% %% Cadenza
% cadenceX=zeros(length(acc),1);
% cadenceY=zeros(length(acc),1);
% tx1=1;
% tx2=0;
% ty1=1;
% ty2=0;
% 
% for i=2:fine
%     if(filteredAcc(i-1,1)>=0 && filteredAcc(i,1)<0)
%         tx2=i;
%         cadenceX(i)=1/((tx2-tx1)*0.04);
%         tx1=i;
%         % disp("rotazione pedale completata in "+num2str(1/cadenceX(i))+"s");
%     end
% 
%     if(filteredAcc(i-1,2)>=0 && filteredAcc(i,2)<0)
%         ty2=i;
%         cadenceY(i)=1/((ty2-ty1)*0.04);
%         ty1=i;
%         % disp("rotazione pedale completata in "+num2str(1/cadenceY(i))+"s");
%     end
% end
% 
% for i=fine-1:-1:1
%     if(cadenceX(i)==0)
%         cadenceX(i)=cadenceX(i+1);
%     end
% 
%     % if(cadenceY(i)==0)
%     %     cadenceY(i)=cadenceY(i+1);
%     % end
% end
% 
% figure
% plot(t,cadenceX,LineWidth=1,Color="g");
% title("cadenza")
% hold on
% plot(t,filteredAcc(:,1)*9.81/-gMedio,LineWidth=1,Color="r");
% grid
% xlabel("t(s)");
% ylabel("rotazioni/s - m/s^2");
% legend("cadenza (rotazioni/s)", "accelerazione (m/s^2)")

% %% Variazione della trasformata di fourier nel tempo (intervalli di 10s)
% 
% rilievo=6;
% 
% % import dei dati
% db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;
% 
% % selezione della porzione di dati da estrarre
% inizio=1;
% fine=length(db);
% 
% intervallo=30;
% 
% acc=acc*gzRot;
% acc=highpass(acc,0.01,25);
% vang=vang*gzRot;
% 
% for i=1:ceil(fine/(25*intervallo))
%     start=(i-1)*(25*intervallo)+1;
%     if(i*(25*intervallo)-1>fine)
%         final=fine;
%     else
%         final=i*(25*intervallo);
%     end
% 
%     f = (0:final-start)*25/(final-start);
% 
% 
%     % figure
%     % plot(t(start:final),acc(start:final,1),LineWidth=1,Color="r");
% 
%     facc=fft(acc(start:final,1));
% 
%     figure
%     plot(f,abs(facc),LineWidth=1,Color="r");
%     title("trasformata discreta di fourier accelerazione in X");
%     xlabel("f (Hz)");
%     ylabel("|X''(f)|");
% 
% end

