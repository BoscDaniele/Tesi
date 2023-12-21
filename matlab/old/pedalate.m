clear
close all
clc

%% Rilevamento
path="..\dbdm\pedalate\";
% path="dbdm\nuovoApproccio2\";

% Questi rilievi sono sati fatti in fondo alla via beita con lo scopo di
% controllare se, come mi aspetto, all'aumentare della velocità diminuisce
% il contributo delle pedalate sull'accelerazione della bici, in alcuni
% rilevamenti ho provato a cambiare, suppongo che questo aumenti
% l'accelerazione trasferita alla bicicletta dalle pedalate.

rilievo=5;

%% Sistema di riferimento sensore = sistema di riferimento bicicletta
[gzRot,gMedio] = GZRot(path);


%% Import + impostaioni
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

inizio=1;
fine=length(db);


t=db(inizio:fine,1)*1e-3;
t=t-t(1);

for i=2:length(t)
    intervalloT(i)=t(i)-t(i-1);
end

disp("intervallo minimo di tempo: "+num2str(min(intervalloT(2:end))));
disp("intervallo massimo di tempo: "+num2str(max(intervalloT(2:end))));

acc=db(inizio:fine,2:4)*gzRot;
vang=db(inizio:fine,5:7)*2*pi/360*1e-3;
mag=db(inizio:fine,8:10)*1e-3;

%% filtraggio dati
% acc=lowpass(acc,5,25);
% acc=highpass(acc,0.2,25);

% acc=lowpass(acc,1,25);
% vang=lowpass(vang,1,25);
% mag=lowpass(mag,1,25);

%% Plot accelerazioni
% plotta3(t,acc,"accelerazioni");
% 
% vel=cumsum(acc)*0.04;
% plotta3(t,vel,"veloctià");

movacc=movmean(acc,250);
% plotta3(t,movacc,"media mobile accelerazioni");

newacc=acc-movacc;
% plotta3(t,newacc,"accelerazioni tolta la media mobile");

newvel=cumsum(newacc)*0.04;
% plotta3(t,newvel,"velocità calcolata togliendo la media mobile dalle accelerazioni");

% acc=newacc;

% stampo le accelerazioni ogni 10s
% for j=1:round(fine/250)
%     i=j-1;
%     if (i*250+250)>fine
%         plotta3(t(round(i*250)+1:fine),acc(round(i*250)+1:fine,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(end)));
%     else
%         plotta3(t(round(i*250)+1:round(i*250)+251),acc(round(i*250)+1:round(i*250)+251,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(round(i*250)+251)));
%     end
% end

%% Plot velocità angolare
% plotta3(t,vang,"velocità angolari");

movevang=movmean(vang,250);
newvang=vang-movevang;

% plotta3(t,newvang,"velocità angolari tolta la media mobile");

% % stampo le accelerazioni ogni 10s
% for j=1:round(fine/250)
%     i=j-1;
%     if (i*250+250)>fine
%         plotta3(t(round(i*250)+1:fine),vang(round(i*250)+1:fine,:),"velocità angolari dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(end)));
%     else
%         plotta3(t(round(i*250)+1:round(i*250)+251),vang(round(i*250)+1:round(i*250)+251,:),"velocità angolari dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(round(i*250)+251)));
%     end
% end


% %% Calcolo velocità, posizione e rotazione mediante integrale
% 
% vel=cumsum(acc)*0.04;
% pos=cumsum(vel)*0.04;
% ang=cumsum(vang)*0.04;
% 
% plotta3(t,ang,"Angoli");
% 
% 
% multiPlotta3(t,acc,vel,"accelerazione","velocità");
% multiPlotta3(t,vel,pos, "velocità","posizione");
% multiPlotta3(t,vang,ang,"velocità angolare", "angoli");

% %% Aggiusto l'orientamento della bicicletta con gli angoli calcolati
%
% for i=1:length(acc)
%     mRot=RotMat(ang(i,:));
%     accr(i,:)=acc(i,:)*mRot;
% end
%
% velR=cumsum(accr)*0.04;
% posR=cumsum(velR)*0.04;
%
% multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
% multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
% multiPlotta3(t,pos,posR, "posizione","posizione ruotata");

%% Prova fourier

% acc=highpass(acc,0.2,25);
% acc=lowpass(acc,25,25);

facc=fft(acc(:,1));
f = (0:length(acc)-1)*25/length(acc);
% f=25/length(acc)*(0:(length(acc)/2));

% p2=abs(facc/length(acc));
% p1=p2(1:length(acc)/2+1);
% p1(2:end-1)=2*p1(2:end-1);

macc=abs(facc);
% macc=abs(p1);

sacc=macc.^2/length(acc);

figure
plot(f,facc,LineWidth=1,Color="r");
title("trasformata discreta di fourier accelerazione in X");
xlabel("Hz");
ylabel("Acc(f)");

figure
plot(f,sacc,LineWidth=1,Color="r");
title("spettro di potenza accelerazione in X");
xlabel("Hz");
ylabel("Sacc(f)");


% froll=fft(vang(:,1));
% mroll=abs(froll);
% % macc=abs(p1);
% 
% sroll=mroll.^2/length(vang);
% 
% figure
% plot(t,vang(:,1),LineWidth=1,Color="r");
% title("velocità angolo di rollio");
% xlabel("s");
% ylabel("rad/s");
% 
% figure
% plot(f,froll,LineWidth=1,Color="r");
% title("trasformata discreta di fourier velocità angolo di rollio");
% xlabel("Hz");
% ylabel("Acc(f)");
% 
% figure
% plot(f,sroll,LineWidth=1,Color="r");
% title("spettro di potenza velocità angolo di rollio");
% xlabel("Hz");
% ylabel("Sacc(f)");

