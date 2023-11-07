close
clear
clc

dbA=importdata("dbdm\Acceleration.csv");
acc=dbA.data;

t=str2num(cell2mat(dbA.textdata(4:end,4)));
t=t-t(1);

%xy
[a,alpha]=CalcolaAngolo(acc(1,1),acc(1,2));

%yz
[b,beta]=CalcolaAngolo(acc(1,2),acc(1,3));

%xz
[c,gamma]=CalcolaAngolo(acc(1,1),acc(1,3));

newX=a*cos(alpha)-c*cos(gamma);
newY=a*sin(alpha)-b*cos(beta);
newZ=b*sin(beta)-c*sin(gamma);
% sqrt(acc(1,1)^2+acc(1,2)^2+acc(1,3)^2);


% figure
% subplot(3,2,1)
% plot([0,acc(1,1)],[0,acc(1,2)],LineWidth=1,Color="b")
% grid
% title("Originale")
% subtitle("XY",Color="b")
% 
% subplot(3,2,2)
% plot([0,newX],[0,newY],LineWidth=1,Color="b")
% grid
% title("Routato")
% subtitle("newXY",Color="b")
% 
% subplot(3,2,3)
% plot([0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="r")
% grid
% subtitle("YZ",Color="r")
% 
% subplot(3,2,4)
% plot([0,newY],[0,newZ],LineWidth=1,Color="r")
% grid
% subtitle("newYZ",Color="r")
% 
% subplot(3,2,5)
% plot([0,acc(1,1)],[0,acc(1,3)],LineWidth=1,Color="g")
% grid
% subtitle("XZ",Color="g")
% 
% subplot(3,2,6)
% plot([0,newX],[0,newZ],LineWidth=1,Color="g")
% grid
% subtitle("newXZ",Color="g")


% figure
% plot3([0,acc(1,1)],[0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="black");
% hold on
% grid
% plot3([0,acc(1,1)],[0,0],[0,0],LineWidth=1,Color="b");
% plot3([0,0],[0,acc(1,2)],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[0,0],[0,acc(1,3)],LineWidth=1,Color="g");
% 
% figure
% plot3([0,newX],[0,newY],[0,newZ],LineWidth=1,Color="black");
% hold on
% grid
% plot3([0,newX],[0,0],[0,0],LineWidth=1,Color="b");
% plot3([0,0],[0,newY],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[0,0],[0,newZ],LineWidth=1,Color="g");

newAcc = zeros(311,3);

for i=1:length(acc)
    newAcc(i,1)=acc(i,1)*cos(alpha)-acc(i,3)*cos(gamma);
    newAcc(i,2)=acc(i,1)*sin(alpha)-acc(i,2)*cos(beta);
    newAcc(i,3)=acc(i,2)*sin(beta)-acc(i,3)*sin(gamma);
end

subplot(3,2,1);
plot(t,acc(:,1),LineWidth=1,Color="b");
title("Originale");
subtitle("X", Color="b");

subplot(3,2,2);
plot(t,-newAcc(:,1),LineWidth=1,Color="b");
title("MenoOriginale");
subtitle("X", Color="b");

subplot(3,2,3);
plot(t,acc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,4);
plot(t,-newAcc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="r");

subplot(3,2,5);
plot(t,acc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

subplot(3,2,6);
plot(t,-newAcc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="g");

%acc(1,:)=(0,0,g)*matriceDiRotazione