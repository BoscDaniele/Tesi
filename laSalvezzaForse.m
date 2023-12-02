clear
close all
clc

dbg=importdata("dbdm\portico\BlueCoin_Log_N000.csv");
g=dbg.data(2+75:end-75,2:4);
t=dbg.data(2+75:end-75,1)*1e-3;
t=t-t(1);

% figure
% subplot(3,1,1);
% plot(t,g(:,1),LineWidth=1,Color="r")
% grid
% subplot(3,1,2);
% plot(t,g(:,2),LineWidth=1,Color="g")
% grid
% subplot(3,1,3);
% plot(t,g(:,3),LineWidth=1,Color="b")
% grid

gmedio=mean(g);

% figure
% plot3([0,gmedio(:,1)],[0,gmedio(:,2)],[0,gmedio(:,3)],LineWidth=1,Color="black");
% hold on
% grid
% axis equal
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="g");
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="b");

modYZ=sqrt(gmedio(2)^2+gmedio(3)^2);
thetaYZ=-acos(-gmedio(3)/modYZ);
mRotX=[1,0,0;0,cos(thetaYZ),-sin(thetaYZ);0,sin(thetaYZ),cos(thetaYZ)];

newgmedio = gmedio*mRotX;

modXz=sqrt(newgmedio(1)^2+newgmedio(3)^2);
thetaXz=-acos(-newgmedio(3)/modXz);
mRoty=[cos(thetaXz),0,sin(thetaXz);0,1,0;-sin(thetaXz),0,cos(thetaXz)];

mRot=mRotX*mRoty;

% lastgmedio=gmedio*mRot;
% 
% figure
% plot3([0,lastgmedio(:,1)],[0,lastgmedio(:,2)],[0,lastgmedio(:,3)],LineWidth=1,Color="black");
% hold on
% grid
% axis equal
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="g");
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="b");

%%
db=importdata("dbdm\portico\BlueCoin_Log_N004.csv");

inizio=2;
fine=length(db.data);

g=norm(gmedio)*1e-2;

acc=db.data(inizio:fine,2:4)*mRot*g*1e-3;
vang=db.data(inizio:fine,5:7)*2*pi/360*1e-3;
t=db.data(inizio:fine,1)*1e-3;
t=t-t(1);

% acc=lowpass(acc,1,25);
% vang=lowpass(vang,1,25);

vel=cumsum(acc)*0.04;
pos=cumsum(vel)*0.04;
ang=cumsum(vang)*0.04;

for i=1:length(acc)
    ang(i,:)
    RotMat(ang(i,:))
    mRot=RotMat(ang(i,:));
    % acc(i,:)
    % smRot=[mRot(i,:,1);mRot(i,:,2);mRot(i,:,3)]

    %smRot=[mRot(i,:,1);mRot(i,:,2);mRot(i,:,3)]
    accr(i,:)=acc(i,:)*mRot;
end

velR=cumsum(accr)*0.04;
posR=cumsum(velR)*0.04;
% angR=cumsum(vangr)*0.04;

% plotta3(t,acc, "accelerazione")
% plotta3(t,accr, "accelerazione ruotata");
% plotta3(t,ang, "angoli");

plotta3(t,pos,"posizione");

% plotta3(t,velR,"velocità ruotata");
plotta3(t,posR,"posizione ruotata");

