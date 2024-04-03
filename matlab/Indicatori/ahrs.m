clear
close all
clc

%% Import dati

path="..\dati\curvaU_forte\";
[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

%% Impostazioni

% i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% rotazione per far coincidere il sistema di riferimento dell'accelerometro
% con quello della bicicletta
rilievo=4;


%% Estrazione e Correzione Dati
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

inizio=1;
fine=length(db);

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*9.81/-gMedio;
vang_rotta=deg2rad(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1);


t=[t_rotta(1)];
last=1;

acc=[acc_rotta(1,:)];
vang=[vang_rotta(1,:)];
mag=[mag_rotta(1,:)];

for i=2:length(t_rotta)
    % disp(num2str(tP(i)-tempoP(end)))
    if(t_rotta(i)-t(end)>=0.039 && t_rotta(i)-t(end)<0.041)
        t=[t,t_rotta(i)];
        acc=[acc;acc_rotta(i,:)];
        vang=[vang;vang_rotta(i,:)];
        mag=[mag;mag_rotta(i,:)];
        last=i;

        if((i-last)>5)
            disp("lammerda, abbiamo saltato 0.004s "+num2str(i));
            brake;
        end
    end
   
end


diff_t=zeros(length(t),1);
diff_t(1)=t(1);
for i=2:length(t)
    diff_t(i)=t(i)-t(i-1);
end


disp("Piano: max diff_t: "+num2str(max(diff_t(2:end))));
disp("Piano: min diff_t: "+num2str(min(diff_t(2:end))));


% %% Import Dati
% db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;
% 
% % selezione della porzione di dati da estrarre
% inizio=1;
% fine=length(db);
% 
% % estrazione dati tempo e conversione in secondi
% t=db(inizio:fine,1)*1e-3;
% t=t-t(1);
% 
% diff_t=zeros(length(t),1);
% diff_t(1)=t(1);
% for i=2:length(t)
%     diff_t(i)=t(i)-t(i-1);
% end
% 
% disp("max diff_t: "+num2str(max(diff_t(2:end))));
% disp("min diff_t: "+num2str(min(diff_t(2:end))));


%% Accelerazione
% acc=db(inizio:fine,2:4)*9.81/-gMedio;
acc_rot=acc*gzRot;

figure("Name","Accelerazione")
subplot(3,1,1)
plot(t,acc_rot(:,1),LineWidth=1,Color="r")
title("Accelerazione")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,acc_rot(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,acc_rot(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Accelerazione Filtrata
filtered_acc=lowpass(acc_rot,0.5,sr);

figure("Name","Accelerazione Filtrata")
subplot(3,1,1)
plot(t,filtered_acc(:,1),LineWidth=1,Color="r")
title("Accelerazione Filtrata")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% Velocità
vel=cumsum(acc_rot)*0.04;

figure("Name","Velocità")
subplot(3,1,1)
plot(t,vel(:,1),LineWidth=1,Color="r")
title("Velocità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s")
grid


%% Velocità Angolare
% vang=deg2rad(db(inizio:fine,5:7)*1e-3);

figure("Name","Velocità Angolare")
subplot(3,1,1)
plot(t,vang(:,1),LineWidth=1,Color="r")
title("Velocità Angolare")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,vang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,vang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad/s")
grid


%% Integrale Velocità Angolare
ang=cumsum(vang)*0.04;

figure("Name","Integrale Velocità Angolare")
subplot(3,1,1)
plot(t,ang(:,1),LineWidth=1,Color="r")
title("Integrale Velocità Angolare")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad")
grid


%% Campo Magnetico
% Gauss*1e-4 -> Tesla
% mT*1e3 -> µT
% Nel codice del sensore il campo magnetico rispetto all'asse y viene preso
% invertito
% mag=([db(inizio:fine,8),-db(inizio:fine,9),db(inizio:fine,10)]*1e-1);

figure("Name","Campo Magnetico")
subplot(3,1,1)
plot(t,mag(:,1),LineWidth=1,Color="r")
title("Campo Magnetico")
subtitle("X")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,2)
plot(t,mag(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("µT")
grid
subplot(3,1,3)
plot(t,mag(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("µT")
grid


%% ahrs
ahrsFuse=ahrsfilter('SampleRate',sr,'OrientationFormat','quaternion', ...
    'ReferenceFrame','NED');

% ahrsFuse=complementaryFilter('SampleRate',sr,'OrientationFormat','quaternion', ...
%     'ReferenceFrame','NED');

[ahrs_orient,ahrs_vang]=ahrsFuse(acc,vang,mag);


%% G
g=rotateframe(ahrs_orient,[0,0,9.81]);
g_rot=g*gzRot;

figure("Name","Gravità")
subplot(3,1,1)
plot(t,g_rot(:,1),LineWidth=1,Color="r")
title("Gravità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,g_rot(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,g_rot(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


acc_noG=(acc-g)*gzRot;

%% ahrs Accelerazione
filtered_acc_noG=lowpass(acc_noG,0.5,sr);

figure("Name","Accelerazione Filtrata Tolta la Gravità")
subplot(3,1,1)
plot(t,filtered_acc_noG(:,1),LineWidth=1,Color="r")
title("Accelerazione Filtrata Tolta la Gravità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,2)
plot(t,filtered_acc_noG(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s^2")
grid
subplot(3,1,3)
plot(t,filtered_acc_noG(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s^2")
grid


%% ahrs Velocità
vel_noG=cumsum(acc_noG)*0.04;

figure("Name","Velocità Tolta la Gravità")
subplot(3,1,1)
plot(t,vel_noG(:,1),LineWidth=1,Color="r")
title("Velocità Tolta la Gravità")
subtitle("X")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,2)
plot(t,vel_noG(:,2),LineWidth=1,Color="g")
subtitle("Y")
xlabel("t(s)")
ylabel("m/s")
grid
subplot(3,1,3)
plot(t,vel_noG(:,3),LineWidth=1,Color="b")
subtitle("Z")
xlabel("t(s)")
ylabel("m/s")
grid


%% ahrs Velocità Angolare

figure("Name","Velocità Angolare AHRS")
subplot(3,1,1)
plot(t,ahrs_vang(:,1),LineWidth=1,Color="r")
title("Velocità Angolare AHRS")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,2)
plot(t,ahrs_vang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad/s")
grid
subplot(3,1,3)
plot(t,ahrs_vang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad/s")
grid


%% ahrs Integrale Velocità Angolare
ahrs_ang=cumsum(ahrs_vang)*0.04;

figure("Name","Integrale Velocità Angolare AHRS")
subplot(3,1,1)
plot(t,ahrs_ang(:,1),LineWidth=1,Color="r")
title("Integrale Velocità Angolare AHRS")
subtitle("Roll")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,2)
plot(t,ahrs_ang(:,2),LineWidth=1,Color="g")
subtitle("Pitch")
xlabel("t(s)")
ylabel("rad")
grid
subplot(3,1,3)
plot(t,ahrs_ang(:,3),LineWidth=1,Color="b")
subtitle("Yaw")
xlabel("t(s)")
ylabel("rad")
grid

