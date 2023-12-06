clear 
close all
clc

path="dbdm\salita-discesa\";

%% Sistema di riferimento sensore = sistema di riferimento bicicletta
inizioG=1;
fineG=inizioG+675;
rilievo=1;

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;
[gRot,gMedio]=GRot(inizioG,fineG,db);

%% Impostazioni
inizio=1;
fine=length(db);

%% Import e filtraggio dati
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*gRot*9.81/-gMedio;
vang=db(inizio:fine,5:7)*2*pi/360*1e-3;

acc=lowpass(acc,1,25);
vang=lowpass(vang,1,25);

%% Calcolo velocità, posizione e rotazione mediante integrale

vel=cumsum(acc)*0.04;
pos=cumsum(vel)*0.04;
ang=cumsum(vang)*0.04;

multiPlotta3(t,acc,vel,"accelerazione","velocità");
multiPlotta3(t,vel,pos, "velocità","posizione");
multiPlotta3(t,vang,ang,"velocità angolare", "angoli");


%% Aggiusto l'orientamento della bicicletta con gli angoli calcolati

for i=1:length(acc)
    mRot=RotMat(ang(i,:));
    accr(i,:)=acc(i,:)*mRot;
end

velR=cumsum(accr)*0.04;
posR=cumsum(velR)*0.04;

multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
multiPlotta3(t,pos,posR, "posizione","posizione ruotata");


% %% Prova Pietro
% 
% acc=db(inizio:fine,2:4)*gRot*9.81/gMedio;
% vang=db(inizio:fine,5:7)*2*pi/360*1e-3;
% 
% acc=lowpass(acc,1,25);
% vang=lowpass(vang,1,25);
% 
% theta=cumsum(vang(:,1))*0.04;
% phi=atan(acc(:,2)/9.81);
% 
% for i=inizio:fine
%     sup=[1, sin(theta(i))*tan(phi(i)), cos(theta(i))*tan(phi(i)); 0, cos(theta(i)), -sin(theta(i)); 0, sin(theta(i))*sec(phi(i)), cos(theta(i))*sec(phi(i))]*vang(i,:)';
%     vangBici(i,:)=sup';
% end
% 
% angBici=cumsum(vangBici)*0.04;
% multiPlotta3(t,vang,vangBici,"velocità angolare prima", "velocità angolare dopo");
% multiPlotta3(t,ang,angBici,"angoli prima", "angoli dopo");
% 
% 
% 
% 
% 
% for i=1:length(angBici)
%     mRot=RotMat(angBici(i,:));
%     accr(i,:)=acc(i,:)*mRot;
% end
% 
% velR=cumsum(accr)*0.04;
% posR=cumsum(velR)*0.04;
% 
% multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
% multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
% multiPlotta3(t,pos,posR, "posizione","posizione ruotata");

% %% Provo imufilter
% 
% acc=db(inizio:fine,2:4)*9.81/gMedio;
% acc=lowpass(acc,1,25);
% 
% FUSE=ahrsfilter('SampleRate',25,OrientationFormat='Rotation matrix');
% 
% [orientation,angularVelocity] = FUSE(acc,vang);
% % orientation=euler(orientation,'XYZ','frame');
% 
% multiPlotta3(t,vang,angularVelocity, "vang","angularVelocity");
% % multiPlotta3(t, ang, orientation, "ang","orientation");
% 
% for i=1:length(acc)
%     % mRot=RotMat(orientation(i,:));
%     % imuAccr(i,:)=acc(i,:)*mRot;
%     imuAccr(i,:)=acc(i,:)*orientation(:,:,i);
% end
% 
% imuVelR=cumsum(imuAccr)*0.04;
% imuPosR=cumsum(imuVelR)*0.04;
% 
% multiPlotta3(t,acc,imuAccr,"accelerazione", "accelerazione imu ruotata");
% multiPlotta3(t,vel,imuVelR, "velocità","velocità imu ruotata");
% multiPlotta3(t,pos,imuPosR, "posizione","posizione imu ruotata");



