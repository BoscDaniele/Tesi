close all
clear
clc

%% Calcolo della Matrice di Rotazione per uniformare il sistema di riferimento del sensore a quello della bicicletta

dbA=importdata("dbdm\nuovo\0.csv");
acc=dbA.data(250:end-325,2:4);

%estrazione tempi
t=dbA.data(250:end-325,1)*1e-3;
t=t-t(1);

media=mean(acc(1:end,:));

modYZ=sqrt(media(2)^2+media(3)^2);
thetaYZ=-acos(-media(3)/modYZ);
mRotX=[1,0,0;0,cos(thetaYZ),-sin(thetaYZ);0,sin(thetaYZ),cos(thetaYZ)];

newMedia = media*mRotX;

modXz=sqrt(newMedia(1)^2+newMedia(3)^2);
thetaXz=-acos(-newMedia(3)/modXz);
mRoty=[cos(thetaXz),0,sin(thetaXz);0,1,0;-sin(thetaXz),0,cos(thetaXz)];

lastMedia = newMedia*mRoty;

mRot=mRotX*mRoty;

%% tutto il resto
db=importdata("dbdm\nuovo\BlueCoin_Log_N001.csv");

%impostazioni
inizio=2;
fine=length(db.data);
g=sqrt(sum(media.^2))*1e-2;

% proviamo a lasciare le accelerazioni così come sono, usando quindi un
% sistema di riferimento solidale alla bicicletta

acc=db.data(inizio:fine,2:4)*mRot*g*1e-3;

ang=db.data(inizio:fine,5:7)*2*pi/360*mRot*1e-3;

t=db.data(inizio:fine,1)*1e-3;
t=t-t(1);

acc=lowpass(acc,25,1e0);
% acc=highpass(acc,1,1e3);
ang=lowpass(ang,40,1e3);

figure(Name="Accelerazioni");
subplot(3,1,1);
plot(t,acc(:,1),LineWidth=1,Color="r");
grid
xlabel("s");
ylabel("m/s^2")
title("Accelerazioni");
subtitle("X");
subplot(3,1,2);
plot(t,acc(:,2),LineWidth=1,Color="g");
grid
xlabel("s");
ylabel("m/s^2")
subtitle("Y");
subplot(3,1,3);
plot(t,acc(:,3),LineWidth=1,Color="b");
grid
xlabel("s");
ylabel("m/s^2")
subtitle("Z");

vel=cumsum(acc)*0.04;

figure(Name="Velocità");
subplot(3,1,1);
plot(t,vel(:,1),LineWidth=1,Color="r");
grid
xlabel("s");
ylabel("m/s")
title("Velocità");
subtitle("X");
subplot(3,1,2);
plot(t,vel(:,2),LineWidth=1,Color="g");
grid
xlabel("s");
ylabel("m/s")
subtitle("Y");
subplot(3,1,3);
plot(t,vel(:,3),LineWidth=1,Color="b");
grid
xlabel("s");
ylabel("m/s")
subtitle("Z");

theta=asin(-acc(:,1)/g);

for i=1:fine-1
phi1(i)=asin(acc(i,2)/(g*cos(theta(i))));
phi2(i)=acos(acc(i,3)/(g*cos(theta(i))));
end

phi1=phi1';
phi2=phi2';

for i=1:length(phi1)
    componenteg(i,:)=[-g*sin(theta(i)),g*cos(theta(i))*sin(phi1(i)),g*cos(theta(i))*cos(phi1(i))];
end

acc=acc+componenteg;

figure(Name="Accelerazioni");
subplot(3,1,1);
plot(t,acc(:,1),LineWidth=1,Color="r");
grid
xlabel("s");
ylabel("m/s^2")
title("Accelerazioni");
subtitle("X");
subplot(3,1,2);
plot(t,acc(:,2),LineWidth=1,Color="g");
grid
xlabel("s");
ylabel("m/s^2")
subtitle("Y");
subplot(3,1,3);
plot(t,acc(:,3),LineWidth=1,Color="b");
grid
xlabel("s");
ylabel("m/s^2")
subtitle("Z");

vel=cumsum(acc)*0.04;

figure(Name="Velocità");
subplot(3,1,1);
plot(t,vel(:,1),LineWidth=1,Color="r");
grid
xlabel("s");
ylabel("m/s")
title("Velocità");
subtitle("X");
subplot(3,1,2);
plot(t,vel(:,2),LineWidth=1,Color="g");
grid
xlabel("s");
ylabel("m/s")
subtitle("Y");
subplot(3,1,3);
plot(t,vel(:,3),LineWidth=1,Color="b");
grid
xlabel("s");
ylabel("m/s")
subtitle("Z");

%funzione imufilter di matlab

FUSE=imufilter('SampleRate',25,'DecimationFactor',2);
t= (0:2:size(acc,1)-1)/25;

[orientation,angularVelocity] = FUSE(acc,ang);
orientation=eulerd(orientation,'XYZ','frame');


figure(Name="orientation");
subplot(3,1,1);
plot(t,orientation(:,1),LineWidth=1,Color="r");
grid
title("orientation");
subtitle("X");
subplot(3,1,2);
plot(t,orientation(:,2),LineWidth=1,Color="g");
grid
subtitle("Y");
subplot(3,1,3);
plot(t,orientation(:,3),LineWidth=1,Color="b");
grid
subtitle("Z");


% figure
% plot(t,orientation(:,1))
% title('Orientation Estimate')
% xlabel('Time (s)')
% ylabel('Rotation (degrees)')

figure(Name="angularVelocity");
subplot(3,1,1);
plot(t,angularVelocity(:,1),LineWidth=1,Color="r");
grid
title("angularVelocity");
subtitle("X");
subplot(3,1,2);
plot(t,angularVelocity(:,2),LineWidth=1,Color="g");
grid
subtitle("Y");
subplot(3,1,3);
plot(t,angularVelocity(:,3),LineWidth=1,Color="b");
grid
subtitle("Z");




