clear
close all
clc

%% Rilevamento
path="dbdm\nuovoApproccio2\";

% Questi rilievi sono sati fatti in fondo alla via secchia
% rilievo 0: vettore g
% rilievo 1: angolo z
% rilievo 2: secchia a scendere
% da qui i dati vengono salvati una volta ogni quando gli pare al sensore
% rilievo 3: salita
% rilievo 4: discesa + frenata
% rilievo 5: discesa + curva a "zig-zag"
% rilievo 6: andata + ritorno secchia-gromlongo
% rilievo 7: secchia a salire

rilievo=2;

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
acc=lowpass(acc,1,25);
vang=lowpass(vang,1,25);
mag=lowpass(mag,1,25);

%% Plot accelerazioni
% plotta3(t,acc,"accelerazioni");

% movacc=movmean(acc,250);
% plotta3(t,movacc,"media mobile accelerazioni");
% 
% newacc=acc-movacc;
% plotta3(t,newacc,"accelerazioni tolta la media mobile");
% 
% newvel=cumsum(newacc)*0.04;
% plotta3(t,newvel,"velocità calcolata togliendo la media mobile dalle accelerazioni");

% acc=newacc;

% % stampo le accelerazioni ogni 10s
% for j=1:round(fine/250)
%     i=j-1;
%     if (i*250+250)>fine
%         plotta3(t(round(i*250)+1:fine),acc(round(i*250)+1:fine,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(end)));
%     else
%         plotta3(t(round(i*250)+1:round(i*250)+251),acc(round(i*250)+1:round(i*250)+251,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(round(i*250)+251)));
%     end
% end

% %% Plot velocità angolare
% plotta3(t,vang,"velocità angolari");
% 
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
