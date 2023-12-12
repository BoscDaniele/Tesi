clear 
close all
clc

%% Rilevamento
path="dbdm\salita-discesa\";
rilievo=6;


%% Sistema di riferimento sensore = sistema di riferimento bicicletta
inizioG=1;
fineG=inizioG+675;

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;
[gRot,gMedio]=GRot(inizioG,fineG,db);

% g=db(inizioG:fineG,2:4)*gRot;
% gMediov=mean(g);
% 
% figure
% plot3([0,gMediov(:,1)],[0,gMediov(:,2)],[0,gMediov(:,3)],LineWidth=1,Color="black");
% hold on
% grid
% axis equal
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="g");
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="b");



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

% multiPlotta3(t,acc,vel,"accelerazione","velocità");
% multiPlotta3(t,vel,pos, "velocità","posizione");
% multiPlotta3(t,vang,ang,"velocità angolare", "angoli");

vangr=vang*gRot;
angr=cumsum(vangr)*0.04;

% multiPlotta3(t,vang,vangr,"velocità angolare", "velocità angolare ruotata");
multiPlotta3(t,ang,angr,"angoli", "angoli ruotati");


%% Aggiusto l'orientamento della bicicletta con gli angoli calcolati

for i=1:length(acc)
    mRot=RotMat(ang(i,:));
    mRotR=RotMat(angr(i,:));
    accr(i,:)=acc(i,:)*mRot;
    accrR(i,:)=acc(i,:)*mRotR;
end

velR=cumsum(accr)*0.04;
posR=cumsum(velR)*0.04;

velRR=cumsum(accrR)*0.04;
posRR=cumsum(velRR)*0.04;

% multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
% multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
% multiPlotta3(t,pos,posR, "posizione","posizione ruotata");
% multiPlotta3(t,pos,posRR, "posizione","posizione ruotata ruotata");
% multiPlotta3(t,posR,posRR, "posizione ruotata","posizione ruotata ruotata");


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

%% Provo imufilter

mag=db(inizio:fine,2:4)*gRot*1e-3;
mag=lowpass(mag,1,25);

FUSE=imufilter('SampleRate',25,'ReferenceFrame','ENU',OrientationFormat='Rotation matrix');

[orientation,angularVelocity] = FUSE(acc,vangr);
% orientation=euler(orientation,'XYZ','frame');

% multiPlotta3(t,vangr,angularVelocity, "vangr","angularVelocity");
% multiPlotta3(t, angr, orientation, "angr","orientation");

for i=1:length(acc)
    % mRot=RotMat(orientation(i,:));
    % imuAccr(i,:)=acc(i,:)*mRot;

    vangR(i,:)=orientation(:,:,i)*vangr(i,:)';

    gravity=orientation(:,:,i)*[0;0;-9.81];
    imuAccr(i,:)=acc(i,:)-gravity';
end

angR=cumsum(vangR)*0.04;

% multiPlotta3(t,vangr,vangR,"vangr", "vangR");
% multiPlotta3(t,angr,angR,"angr", "angR");

% imuVelR=cumsum(imuAccr)*0.04;
% imuPosR=cumsum(imuVelR)*0.04;
% 
% multiPlotta3(t,acc,imuAccr,"accelerazione", "accelerazione imu ruotata");
% multiPlotta3(t,vel,imuVelR, "velocità","velocità imu ruotata");
% multiPlotta3(t,pos,imuPosR, "posizione","posizione imu ruotata");


for i=1:length(imuAccr)
    mRot=RotMat(-angR(i,:));
    accr(i,:)=imuAccr(i,:)*mRot;
end

velR=cumsum(accr)*0.04;
posR=cumsum(velR)*0.04;


multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
multiPlotta3(t,pos,posR, "posizione","posizione ruotata");
