clear
close all
clc

%% Import dati

path="db\secchia\";
[gzRot,gMedio] = GZRot(path);

% selezionare il rilievo da caricare
% 0 - gravità
% 1 - inclinazione
% 2 - discesa via secchia
% 5 - curve
% 6 - pedalata tranquilla
rilievo=5;

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
% StampaAcc(t,acc,"acc","acc")

% Visto che la gravità non è stata rimossa dalle misure, questa agisce come
% un'accelerazione negativa, quando la bicicletta si sta movendo in discesa,
% o un'accelerazione positiva, quando la bicicletta si sta movendo in
% salita, lungo l'asse delle x. Questo ha senso se si pensa che la gravità
% è un vettore negativo anche se diretto verso il basso (il sensore pensa
% che ci sia una forza pari a g che lo tira costantemente verso l'alto)

% acc_mean=lowpass(acc,1,sr);
% acc_mean=highpass(acc_mean,0.01,sr);
acc_mean=movmean(acc,[40,0]);
% StampaAcc(t,acc_mean,"acc media","acc media")

vang=deg2rad(db(inizio:fine,5:7)*1e-3);
StampaVang(t,vang,"vang","vang")
% ang=cumsum(vang)*0.04;
% StampaAng(t,ang,"ang","ang")


%% PP (pedaling pattern) - ampiezza oscillazione accelerazione lungo x
min_value=0;
max_value=0;

% stato: stato in cui mi trovo nella macchina a stati finiti
% 0: fase crescente dell'accelerazione
% 1: fase calante dell'accelerazione
stato=0;

ampiezza=zeros(length(acc),1);
first_amp=0;

for i=2:length(acc)
    if(~stato && acc(i-1,1)>acc(i,1))
        stato=1;
        max_value=acc(i-1,1);
        ampiezza(i)=max_value-min_value;
    elseif(stato && acc(i-1,1)<acc(i,1))
        stato=0;
        min_value=acc(i-1,1);
        % ampiezza(i)=max_value-min_value;
        if(~first_amp)
            first_amp=ampiezza(i);
        end
    end

    if(ampiezza(i)==0)
        ampiezza(i)=ampiezza(i-1);
    end   
end

figure
plot(t,ampiezza,LineWidth=1,Color="r");
title("Ampiezza Oscillazione");
xlabel("t(s)")
grid
hold on
plot(t,acc(:,1),LineWidth=1,Color="b");
legend("Ampiezza(m/s^2)","accelerazione(m/s^2)")


%% PP (pedaling pattern) - ampiezza oscillazione accelerazione lungo x - accelerazione filtrata
filtered_acc=highpass(acc,1,sr);
filtered_acc=lowpass(filtered_acc,5,sr);

% StampaAcc(t,filtered_acc,"filtered_acc","filtered_acc")

min_value_f=0;
max_value_f=0;

% stato: stato in cui mi trovo nella macchina a stati finiti
% 0: fase crescente dell'accelerazione
% 1: fase calante dell'accelerazione
stato_f=0;

ampiezza_f=zeros(length(acc),1);
first_amp_f=0;

for i=2:length(filtered_acc)
    if(~stato_f && filtered_acc(i-1,1)>filtered_acc(i,1))
        stato_f=1;
        max_value_f=filtered_acc(i-1,1);
        ampiezza_f(i)=max_value_f-min_value_f;

    elseif(stato_f && filtered_acc(i-1,1)<filtered_acc(i,1))
        stato_f=0;
        min_value_f=filtered_acc(i-1,1);
        % ampiezza_f(i)=max_value_f-min_value_f;
        if(~first_amp_f)
            first_amp_f=ampiezza_f(i);
        end
    end

    if(ampiezza_f(i)==0)
        ampiezza_f(i)=ampiezza_f(i-1);
    end
end

figure
plot(t,ampiezza_f,LineWidth=1,Color="r");
title("Ampiezza Oscillazione Accelerazione Filtrata");
xlabel("t(s)")
grid
hold on
plot(t,filtered_acc(:,1),LineWidth=1,Color="b");
legend("Ampiezza(m/s^2)","accelerazione 1-5Hz (m/s^2)")


% %% Confronto ampiezza filtrata e non
% figure
% plot(t,ampiezza,LineWidth=1,Color="r");
% title("Ampiezza Oscillazione Accelerazione Filtrata");
% xlabel("t(s)")
% grid
% hold on
% plot(t,ampiezza_f,LineWidth=1,Color="b");
% legend("Ampiezza(m/s^2)","Ampiezza 1-5Hz (m/s^2)")


%% Confronto Accelerazione Media Ampiezza
zero=zeros(length(acc),1);

figure
plot(t,ampiezza,LineWidth=1,Color="r")
title("Accelerazione Media - Ampiezza")
xlabel("t(s)")
ylabel("accelerazione(m/s^2)")
grid
hold on
plot(t,acc_mean(:,1),LineWidth=1,Color="b")
plot(t,zero,Color="black")
legend("Ampiezza","Accelerazione Media","0")

figure
plot(t,ampiezza_f,LineWidth=1,Color="r")
title("Accelerazione Media - Ampiezza Filtrata")
xlabel("t(s)")
ylabel("accelerazione(m/s^2)")
grid
hold on
plot(t,acc_mean(:,1),LineWidth=1,Color="b")
plot(t,zero,Color="black")
legend("Ampiezza Filtrata","Accelerazione Media","0")


% %% EPF
% disp("ampiezza prima oscillazione: "+num2str(first_amp));
% disp("ampiezza prima oscillazione: "+num2str(first_amp_f)+" (accelerazione filtrata)");
% 
% figure
% plot(t,ampiezza/first_amp,LineWidth=1,Color="r");
% title("Ampiezza Oscillazione/Prima Ampiezza");
% xlabel("t(s)")
% ylabel("ampiezza/prima ampiezza")
% grid
% 
% figure
% plot(t,ampiezza_f/first_amp_f,LineWidth=1,Color="r");
% title("Ampiezza Oscillazione/Prima Ampiezza Filtrata");
% xlabel("t(s)")
% ylabel("ampiezza/prima ampiezza")
% grid


%% Cadenza accelerazione
stato=0;
rot_period=zeros(length(acc),1);
t_min=0;
t_max=0;

for i=2:length(acc)
    if(~stato && acc(i,1)<acc(i-1,1))
        t_max=i;
        rot_period(i)=(t_max-t_min)*0.04;
        stato=1;
    elseif(stato && acc(i,1)>acc(i-1,1))
        t_min=i;
        stato=0;
    end

    if(rot_period(i)==0)
        rot_period(i)=rot_period(i-1);
    end
end

cadence=rot_period.^-1;

figure
plot(t,cadence,LineWidth=1,Color="r")
title("Cadenza")
xlabel("t(s)")
grid
hold on
plot(t,acc(:,1),LineWidth=1,Color="b")
legend("cadence(1/s)","accelerazione")


%% Cadenza accelerazione filtrata
stato_f=0;
rot_period_f=zeros(length(filtered_acc),1);
t_min_f=0;
t_max_f=0;

for i=2:length(filtered_acc)
    if(~stato_f && filtered_acc(i,1)<filtered_acc(i-1,1))
        t_max_f=i;
        rot_period_f(i)=(t_max_f-t_min_f)*0.04;
        stato_f=1;
    elseif(stato_f && filtered_acc(i,1)>filtered_acc(i-1,1))
        t_min_f=i;
        stato_f=0;
    end

    if(rot_period_f(i)==0)
        rot_period_f(i)=rot_period_f(i-1);
    end
end

cadence_f=rot_period_f.^-1;

figure
plot(t,cadence_f,LineWidth=1,Color="r")
title("Cadenza Accelerazione Filtrata")
xlabel("t(s)");
grid
hold on
plot(t,filtered_acc(:,1),LineWidth=1,Color="b")
legend("cadence(1/s)","accelerazione")


%% Confronto Cadenza Filtratata e non
figure
plot(t,cadence,LineWidth=1,Color="r")
title("Confronto Cadenze")
xlabel("t(s)");
grid
hold on
plot(t,cadence_f,LineWidth=1,Color="b")
legend("Cadenza", "Cadenza Filtrata")


%% Confronto Ampiezza-Cadenza
figure
plot(t,ampiezza,LineWidth=1,Color="b")
title("Ampiezza-Cadenza")
xlabel("t(s)")
grid
hold on
plot(t,cadence,LineWidth=1,Color="r")
legend("Ampiezza","Cadenza")

figure
plot(t,ampiezza_f,LineWidth=1,Color="b")
title("Ampiezza-Cadenza Filtrate")
xlabel("t(s)")
grid
hold on
plot(t,cadence_f,LineWidth=1,Color="r")
legend("Ampiezza","Cadenza")


% Affinchè funzioni la rimozione della gravità è necessario capire come
% ottenere la rotazione che mi riporta le misure nel sistema di riferimento
% originario, perchè, ovviamente, invertendo la rotazione che la funzione
% ahrsfilter mi restituisce non ottengo la rotazione da applicare per
% tornare indietro

% %% prova ahrsfilter per eliminare g
% acc=db(inizio:fine,2:4)*9.81/-gMedio;
% StampaAcc(t,acc,"acc","acc no rot")
% 
% mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);
% 
% fuse = ahrsfilter('SampleRate',sr,'OrientationFormat','Rotation matrix', ...
%     'ReferenceFrame','NED');
% 
% [orientationMatrix,~] = fuse(acc,vang,mag);
% 
% orientationMatrix_inv=zeros(3,3,length(orientationMatrix));
% for i=1:length(orientationMatrix_inv)
%     orientationMatrix_inv(:,:,i)=orientationMatrix(:,:,i)';
% end
% 
% % % ruota g
% % g=zeros(length(acc),3);
% % for i=1:length(acc)
% %     g(i,:)=[0,0,9.81]*orientationMatrix_inv(:,:,i);
% % end
% % StampaAcc(t,g,"g","g")
% % 
% % newAcc=acc-g;
% % StampaAcc(t,newAcc,"newAcc","newAcc no rot")
% % 
% % StampaAcc(t,acc*gzRot,"acc","acc rot")
% % StampaAcc(t,newAcc*gzRot,"newAcc rot","newAcc rot")
% % 
% % acc_medio=movmean(acc*gzRot,[40,0]);
% % newAcc_medio=movmean(newAcc*gzRot,[40,0]);
% % 
% % StampaAcc(t,acc_medio,"acc medio","acc medio")
% % StampaAcc(t,newAcc_medio,"newAcc medio","newAcc medio")
% 
% % % ruota acc
% % newAcc=zeros(length(acc),3);
% % for i=1:length(acc)
% %     newAcc(i,:)=acc(i,:)*orientationMatrix(:,:,i);
% % end
% % StampaAcc(t,newAcc,"newAcc","newAcc")
% % 
% % newAcc=newAcc-[0,0,9.81];
% % StampaAcc(t,newAcc,"newAcc no g","newAcc no g")
% % 
% % accNoG=zeros(length(acc),3);
% % for i=1:length(acc)
% %     accNoG(i,:)=newAcc(i,:)*orientationMatrix_inv(:,:,i);
% % end
% % StampaAcc(t,accNoG,"accNoG","accNoG")
% % 
% % StampaAcc(t,acc*gzRot,"acc","acc rot")
% % StampaAcc(t,accNoG*gzRot,"accNoG rot","accNoG rot")
% % 
% % acc_medio=movmean(acc*gzRot,[40,0]);
% % accNoG_medio=movmean(accNoG*gzRot,[40,0]);
% % 
% % StampaAcc(t,acc_medio,"acc medio","acc medio")
% % StampaAcc(t,accNoG_medio,"accNoG medio","accNoG medio")


% vel=cumsum(acc)*0.04;
% StampaAcc(t,acc_mean,"acc mean","acc mean")
% StampaVel(t,vel,"vel","vel")


%% prova norma
norma=zeros(length(acc),1);
for i=1:length(acc)
    norma(i)=norm(acc(i,:))-9.81;
end

norma_media=movmean(norma,[40,0]);

figure
plot(t,norma,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,norma_media,LineWidth=1,Color="b")
legend("vettore accelerazione(m/s^2)","media vettore accelerazione(m/s^2)")


figure
plot(t,ampiezza_f,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,norma_media,LineWidth=1,Color="b")
legend("ampiezza(m/s^2)","norma vettore accelerazione(m/s^2)")


normaXY=zeros(length(acc),1);
normaXZ=zeros(length(acc),1);
normaYZ=zeros(length(acc),1);

for i=1:length(acc)
    normaXY(i)=norm(acc(i,1:2));
    normaXZ(i)=norm([acc(i,1),acc(i,3)])-9.81;
    normaYZ(i)=norm(acc(i,2:3))-9.81;
end

normaXY_medio=movmean(normaXY,[40,0]);
normaXZ_medio=movmean(normaXZ,[40,0]);
normaYZ_medio=movmean(normaYZ,[40,0]);

figure
plot(t,normaXY_medio,LineWidth=1,Color="r")
xlabel("t")
grid
hold on
plot(t,normaXZ_medio,LineWidth=1,Color="g")
plot(t,normaYZ_medio,LineWidth=1,Color="b")
legend("normaXY","normaXZ","normaYZ")


figure
plot(t,ampiezza_f,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,normaXY_medio,LineWidth=1,Color="b")
legend("ampiezza(m/s^2)","media vettore accelerazione(m/s^2)")

figure
plot(t,filtered_acc,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,normaXY_medio,LineWidth=1,Color="b")
plot(t,vang(:,3),LineWidth=1,Color="black")
legend("accelerazione filtrata(m/s^2)","media vettore accelerazione(m/s^2)","velocità angolare Z(rad/s)")

figure
plot(t,filtered_acc,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,vang(:,1),LineWidth=1,Color="black")
legend("accelerazione filtrata(m/s^2)","velocità angolare X(rad/s)")

figure
plot(t,filtered_acc,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,vang(:,2),LineWidth=1,Color="black")
legend("accelerazione filtrata(m/s^2)","velocità angolare Y(rad/s)")

figure
plot(t,filtered_acc,LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,vang(:,3),LineWidth=1,Color="black")
legend("accelerazione filtrata(m/s^2)","velocità angolare Z(rad/s)")


figure
plot(t,vang(:,3),LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,norma_media,LineWidth=1,Color="b")
legend("velocità angolare(rad/s)","media vettore accelerazione(m/s^2)")


figure
plot(t,vang(:,3),LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,normaXY_medio,LineWidth=1,Color="b")
legend("velocità angolare(rad/s)","norma vettore accelerazione XY(m/s^2)")


figure
plot(t,vang(:,3),LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,normaXZ_medio,LineWidth=1,Color="b")
legend("velocità angolare(rad/s)","norma vettore accelerazione XZ(m/s^2)")


figure
plot(t,vang(:,3),LineWidth=1,Color="r")
xlabel("t(s)")
grid
hold on
plot(t,normaYZ_medio,LineWidth=1,Color="b")
legend("velocità angolare(rad/s)","norma vettore accelerazione YZ(m/s^2)")


%% Prova Varianza
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
plot(t,v_acc(:,1),LineWidth=1,Color="black")
xlabel("t(s)")
grid
hold on
plot(t,ampiezza,LineWidth=1,Color="r")
legend("varianza","ampiezza")




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
