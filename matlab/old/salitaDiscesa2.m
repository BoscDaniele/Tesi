clear
close all
clc

%% Rilevamento
path="dbdm\salita-discesa2\";
rilievo=2; % il rilievo 0 è la gravità, il rilievo 1 è per stabilire di quanto è ruotato il sensore rispetto all'asse z


%% Sistema di riferimento sensore = sistema di riferimento bicicletta
[gzRot,gMedio] = GZRot(path);

%% Import + impostaioni
db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

inizio=1;
fine=length(db);


t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*gzRot*9.81/-gMedio;
vang=db(inizio:fine,5:7)*gzRot*2*pi/360*1e-3;
mag=db(inizio:fine,8:10)*gzRot*1e-3;


%% filtraggio dati
acc=lowpass(acc,1,25);
vang=lowpass(vang,1,25);
mag=lowpass(mag,1,25);

%% Calcolo velocità, posizione e rotazione mediante integrale

vel=cumsum(acc)*0.04;
pos=cumsum(vel)*0.04;
ang=cumsum(vang)*0.04;

% multiPlotta3(t,acc,vel,"accelerazione","velocità");
% multiPlotta3(t,vel,pos, "velocità","posizione");
% multiPlotta3(t,vang,ang,"velocità angolare", "angoli");

%% Aggiusto l'orientamento della bicicletta con gli angoli calcolati

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

%% Provo imufilter
FUSE=imufilter('SampleRate',25,'ReferenceFrame','ENU',OrientationFormat='Rotation matrix');

[orientation,angularVelocity] = FUSE(acc,vang);
% orientation=euler(orientation,'XYZ','frame');

% multiPlotta3(t,vang,angularVelocity, "vang","angularVelocity");
% multiPlotta3(t, ang, orientation, "ang","orientation");

for i=1:length(acc)
    % mRot=RotMat(orientation(i,:));
    % imuAccr(i,:)=acc(i,:)*mRot;

    % imuAccr(i,:)=acc(i,:)*orientation(:,:,i);
    vangR(i,:)=vang(i,:)*orientation(:,:,i);

    gravity=orientation(:,:,i)*[0;0;-9.81];
    imuAccr(i,:)=(orientation(:,:,i)*acc(i,:)'-gravity)';

    % gravity=[0,0,-9.81]*orientation(:,:,i);
    % imuAccr(i,:)=acc(i,:)*orientation(:,:,i)-gravity;
end

angR=cumsum(vangR)*0.04;

% multiPlotta3(t,vang,vangR,"vang", "vangR");
% multiPlotta3(t,ang,angR,"ang", "angR");

imuVelR=cumsum(imuAccr)*0.04;
imuPosR=cumsum(imuVelR)*0.04;

multiPlotta3(t,acc,imuAccr,"accelerazione", "accelerazione imu ruotata");
multiPlotta3(t,vel,imuVelR, "velocità","velocità imu ruotata");
multiPlotta3(t,pos,imuPosR, "posizione","posizione imu ruotata");


% for i=1:length(imuAccr)
%     mRot=RotMat(angR(i,:));
%     accr(i,:)=imuAccr(i,:)*mRot;
% end
%
% velR=cumsum(accr)*0.04;
% posR=cumsum(velR)*0.04;
%
%
% multiPlotta3(t,acc,accr,"accelerazione","accelerazione ruotata");
% multiPlotta3(t,vel,velR, "velocità","velocità ruotata");
% multiPlotta3(t,pos,posR, "posizione","posizione ruotata");
%



