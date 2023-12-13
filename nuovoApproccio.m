clear
close all
clc

% Questo nuovo approccio parte male... Le misure per resettare la gravità e
% la rotazione attorno all'asse z non ci sono

%% Rilevamento
path="dbdm\nuovoApproccio\";

% Questi rilievi sono sati fatti in fondo alla via secchia
% rilievo 0: discesa
% rilievo 1: salita
% rilievo 2: discesa
% rilievo 3: dalla via secchia alla via beita
% rilievo 4: dalla via beita alla via secchia
% rilievo 5: salita
% rilievo 6: parte "piana" della via secchia

rilievo=1;

%% Import + impostaioni
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

inizio=1;
fine=length(db);


t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4);
vang=db(inizio:fine,5:7)*2*pi/360*1e-3;
mag=db(inizio:fine,8:10)*1e-3;

%% Plot accelerazioni
plotta3(t,acc,"accelerazioni");

% stampo le accelerazioni ogni 10s
for j=1:round(fine/250)
    i=j-1;
    if (i*250+250)>fine
        plotta3(t(round(i*250)+1:fine),acc(round(i*250)+1:fine,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(end)));
    else
        plotta3(t(round(i*250)+1:round(i*250)+251),acc(round(i*250)+1:round(i*250)+251,:),"accelerazioni dal "+num2str(t(round(i*250)+1))+" secondo al "+num2str(t(round(i*250)+251)));
    end
end


%% filtraggio dati
acc=lowpass(acc,1,25);
vang=lowpass(vang,1,25);
mag=lowpass(mag,1,25);

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

