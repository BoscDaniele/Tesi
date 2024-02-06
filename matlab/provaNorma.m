clear
close all
clc

%% Import dati

path="db\palazzago2\";
[gzRot,gMedio] = GZRot(path);

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - discesa + curve + alla fine salita
% 3 - frenata posteriore
% 4 - frenata anteriore
% 5 - percorso misto (leggera discesa + rotonda + leggera salita)
% 6 - discesa + molta salita
rilievo=7;

sr = 25; %sample rate


%% Import Dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

% selezione della porzione di dati da estrarre
inizio=1;
fine=length(db);

% estrazione dati tempo e conversione in secondi
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*gzRot*9.81/-gMedio;
StampaAcc(t,acc,"acc","acc")

acc_mean=movmean(acc,[40,0]);
StampaAcc(t,acc_mean,"acc media","acc media")

vang=deg2rad(db(inizio:fine,5:7)*1e-3);
StampaVang(t,vang,"vang","vang")


%% Norma Vettore Accelerazione
norma_Acc=zeros(length(acc),1);
for i=1:length(acc)
    % norma_Acc(i)=norm(acc(i,:));
    norma_Acc(i)=norm([acc(i,1),acc(i,2)]);
end

figure
plot(t,norma_Acc,LineWidth=1,Color="r")
title("Norma Accelerazione")
xlabel("t(s)")
ylabel("Accelerazione(m/s^2)")
grid


mean_norma_Acc=movmean(norma_Acc,[40,0]);

figure
plot(t,mean_norma_Acc,LineWidth=1,Color="r")
title("Norma Accelerazione Media")
xlabel("t(s)")
ylabel("Accelerazione(m/s^2)")
grid


%% ahrsfilter

acc=db(inizio:fine,2:4)*9.81/-gMedio;
vang=deg2rad(db(inizio:fine,5:7)*1e-3);
% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);

fuse=ahrsfilter('SampleRate',sr,'OrientationFormat','Rotation matrix', ...
    'ReferenceFrame','NED');
[orientation,angularVelocity] = fuse(acc,vang,mag);

newG=zeros(length(acc),3);
for i=1:length(newG)
    newG(i,:)=[0,0,9.81]*orientation(:,:,i)';
end

StampaAcc(t,newG,"newG","newG")

newAcc=(acc-newG)*gzRot;
StampaAcc(t,newAcc,"newAcc","newAcc")

newAcc_mean=movmean(newAcc,[40,0]);
StampaAcc(t,newAcc_mean,"media newAcc","media newAcc")


norma_newAcc=zeros(length(acc),1);
for i=1:length(acc)
    % norma_Acc(i)=norm(acc(i,:));
    norma_newAcc(i)=norm([newAcc(i,:)]);
end

figure
plot(t,norma_newAcc,LineWidth=1,Color="r")
title("Norma Accelerazione")
xlabel("t(s)")
ylabel("Accelerazione(m/s^2)")
grid


mean_norma_newAcc=movmean(norma_newAcc,[40,0]);

figure
plot(t,mean_norma_newAcc,LineWidth=1,Color="r")
title("Norma Accelerazione Media")
xlabel("t(s)")
ylabel("Accelerazione(m/s^2)")
grid


v_acc=movvar(acc,[40,0]);

figure
plot(t,v_acc(:,1),LineWidth=1,Color="r")
title("Varianza")
xlabel("t(s)")
ylabel("varianza")
grid
hold on
plot(t,v_acc(:,2),LineWidth=1,Color="g")
plot(t,v_acc(:,3),LineWidth=1,Color="b")
legend("varianza X","varianza Y","varianza Z")


figure
plot(t,v_acc(:,1),LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,mean_norma_Acc,LineWidth=1,Color="b")
legend("varianza","norma")



v_norma=movvar(norma_Acc,[40,0]);

figure
plot(t,v_norma,LineWidth=1,Color="r");
title("Varianza Norma")
xlabel("t(s)")
ylabel("varianza")
grid

figure
plot(t,v_acc(:,1),LineWidth=1,Color="r")
title("Varianza")
xlabel("t(s)")
ylabel("varianza")
grid
hold on
plot(t,v_acc(:,2),LineWidth=1,Color="g")
plot(t,v_acc(:,3),LineWidth=1,Color="b")
plot(t,v_norma,LineWidth=1,Color="black");
legend("varianza X","varianza Y","varianza Z","varianza Norma")






%% Funzioni
function Stampa(x,y,nome,titolo,sottotitolo,etichettaX,etichettaY)
figure(Name=nome)
subplot(3,1,1)
plot(x,y(:,1),LineWidth=1,Color="r");
title(titolo)
subtitle(sottotitolo(1))
grid
xlabel(etichettaX);
ylabel(etichettaY(1));

subplot(3,1,2)
plot(x,y(:,2),LineWidth=1,Color="g");
subtitle(sottotitolo(2))
grid
xlabel(etichettaX);
ylabel(etichettaY(2));

subplot(3,1,3)
plot(x,y(:,3),LineWidth=1,Color="b");
subtitle(sottotitolo(3))
grid
xlabel(etichettaX);
ylabel(etichettaY(3));
end

function StampaAcc(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X''(m/s^2)","Y''(m/s^2)","Z''(m/s^2)"])
end

function StampaVel(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X'(m/s)","Y'(m/s)","Z'(m/s)"])
end

function StampaPos(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"t(s)",["X(m)","Y(m)","Z(m)"])
end

function StampaVang(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"t(s)",["r'(rad/s)","p'(rad/s)","y'(rad/s)"])
end

function StampaAng(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"t(s)",["r(rad)","p(rad)","y(rad)"])
end

function StampaFreqAcc(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["X","Y","Z"],"f(Hz)",["|X''(f)|","|Y''(f)|","|Z''(f)|"])
end

function StampaFreqVang(x,y,nome,titolo)
Stampa(x,y,nome,titolo,["Roll","Pitch","Yaw"],"f(Hz)",["|r'(f)|","|p'(f)|","|y'(f)|"])
end

