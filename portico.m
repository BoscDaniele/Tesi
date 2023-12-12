clear
close all
clc

dbg=importdata("dbdm\portico\BlueCoin_Log_N001.csv");
g=dbg.data(2+150:end-75,2:4);
t=dbg.data(2+150:end-75,1)*1e-3;
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
db1=importdata("dbdm\portico\BlueCoin_Log_N007.csv");

inizio=2;
fine=length(db1.data);

g=norm(gmedio)*1e-2;

acc=db1.data(inizio:fine,2:4)*g*1e-3;
ang=db1.data(inizio:fine,5:7)*mRot*(2*pi/360)*1e-3;
t=db1.data(inizio:fine,1)*1e-3;
t=t-t(1);

% figure
% subplot(3,1,1);
% plot(t,acc(:,1),LineWidth=1,Color="r")
% grid
% subplot(3,1,2);
% plot(t,acc(:,2),LineWidth=1,Color="g")
% grid
% subplot(3,1,3);
% plot(t,acc(:,3),LineWidth=1,Color="b")
% grid

% acc=lowpass(acc,25,1e0);
accr=acc*mRot;

figure
subplot(3,1,1);
plot(t,ang(:,1),LineWidth=1,Color="r")
grid
subplot(3,1,2);
plot(t,ang(:,2),LineWidth=1,Color="g")
grid
subplot(3,1,3);
plot(t,ang(:,3),LineWidth=1,Color="b")
grid

vel=cumsum(accr)*0.04;

figure(Name="velocità")
subplot(3,1,1);
plot(t,vel(:,1),LineWidth=1,Color="r")
grid
subplot(3,1,2);
plot(t,vel(:,2),LineWidth=1,Color="g")
grid
subplot(3,1,3);
plot(t,vel(:,3),LineWidth=1,Color="b")
grid

spaz=cumsum(vel)*0.04;

figure(Name="spazio");
subplot(3,1,1);
plot(t,spaz(:,1),LineWidth=1,Color="r")
grid
subplot(3,1,2);
plot(t,spaz(:,2),LineWidth=1,Color="g")
grid
subplot(3,1,3);
plot(t,spaz(:,3),LineWidth=1,Color="b")
grid

% accM=accr-mean(accr);
% 
% figure
% subplot(3,1,1);
% plot(t,accM(:,1),LineWidth=1,Color="r")
% grid
% subplot(3,1,2);
% plot(t,accM(:,2),LineWidth=1,Color="g")
% grid
% subplot(3,1,3);
% plot(t,accM(:,3),LineWidth=1,Color="b")
% grid
% 
% velM=cumsum(accM)*0.04;
% 
% figure(Name="velocità")
% subplot(3,1,1);
% plot(t,velM(:,1),LineWidth=1,Color="r")
% grid
% subplot(3,1,2);
% plot(t,velM(:,2),LineWidth=1,Color="g")
% grid
% subplot(3,1,3);
% plot(t,velM(:,3),LineWidth=1,Color="b")
% grid
% 
% spazM=cumsum(velM)*0.04;
% 
% figure(Name="spazio");
% subplot(3,1,1);
% plot(t,spazM(:,1),LineWidth=1,Color="r")
% grid
% subplot(3,1,2);
% plot(t,spazM(:,2),LineWidth=1,Color="g")
% grid
% subplot(3,1,3);
% plot(t,spazM(:,3),LineWidth=1,Color="b")
% grid

