clear all
close all
clc

dati={{".\dati\lunga_forte\",2},{".\dati\curvaU_forte\",4},{".\dati\curva_forte\",4},{".\db\secchia\",2}};

n=1;

path=dati{n}{1};
rilievo=dati{n}{2};

% path=".\dati\curvaU_piano\";
% path=".\db\secchia\";
% rilievo=2;

sr=25; % sample rate

[gzRot,gMedio] = GZRot(path);

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*gzRot*9.81/-gMedio;
vang_rotta=(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1)*gzRot;

[t,acc,vang,mag]=AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta);
vel=cumsum(acc)*0.04;

font="Times New Roman";

% Rettilineo
linee=ones(5,1);
linee(1)=1;
linee(2)=15;
linee(3)=19.4;
linee(4)=23.76;
linee(5)=25.4;

str={'Accelerate','Idle','Accelerate','Idle','Brake + Turning'};


% % Curva U
% linee=ones(4,1);
% linee(1)=2;
% linee(2)=15;
% linee(3)=20;
% linee(4)=28;
% 
% str={'Accelerate 1','Turning','Accelerate 2','Brake'};

fun=acc;

axs={"X","Y","Z"};
xlbl="t(s)";
ylbl="Acc. (m/s^2)";

% axs={"Roll","Pitch","Yaw"};
% xlbl="t(s)";
% ylbl="V. Ang. (deg/s)";

% axs={"X","Y","Z"};
% xlbl="t(s)";
% ylbl="Mag. Field (µT)";



limY=zeros(3,2);
limY(1,:)=[floor(min(fun(:,1))),ceil(max(fun(:,1)))];
limY(2,:)=[floor(min(fun(:,2))),ceil(max(fun(:,2)))];
limY(3,:)=[floor(min(fun(:,3))),ceil(max(fun(:,3)))];

limY(2,:)=[min([limY(1,:),limY(2,:)]),max([limY(1,:),limY(2,:)])];


textLimY=ones(3,5);

for i=1:3
    textLimY(i,:)=(limY(i,2)+abs(limY(i,2)).*0.1)*ones(5,1);
end

limY(:,2)=(limY(:,2)+1)+abs(limY(:,2)+1).*0.25;
% limY(3,2)=limY(3,2)+1;

% textLimY=ones(3,4);

% for i=1:3
%     textLimY(i,:)=(limY(i,2)+abs(limY(i,2)).*0.1)*ones(4,1);
% end
% 
% limY(:,2)=(limY(:,2)+1)+abs(limY(:,2)+1).*0.25;
% limY(3,2)=limY(3,2)+1;

% textLimY(1,:)=(limY(1,2)+abs(limY(1,2)).*0.1)*ones(4,1);
% textLimY(2,:)=(limY(2,2)+abs(limY(2,2)).*0.1)*ones(4,1);
% textLimY(3,:)=(40)*ones(4,1);
% 
% limY(1,2)=(limY(1,2)+1)+abs(limY(1,2)+1).*0.25;
% limY(2,2)=(limY(2,2)+1)+abs(limY(2,2)+1).*0.25;
% limY(3,2)=(limY(3,2)+10)+abs(limY(3,2)+10).*0.25;



f=figure;
subplot(3,1,1)
plot(t,fun(:,1),LineWidth=1,color="r")
title("Acceleration",FontName=font)
subtitle(axs(1),FontName=font)
xlabel(xlbl,FontName=font)
ylabel(ylbl,FontName=font)
ylim(limY(1,:))
grid
hold on
for j=1:length(linee)
    line(linee(j)*[1,1],limY(1,:), 'Color','black')
end
text(linee(:)+0.5*[1,1,1,1/2,1]',textLimY(1,:),str,FontSize=6.5, FontWeight="bold",FontName=font)
% text(linee(:)+0.5*[1,1,1,1]',textLimY(1,:),str,FontSize=6.5, FontWeight="bold",FontName=font)

subplot(3,1,2)
plot(t,fun(:,2),LineWidth=1,color="g")
subtitle(axs(2),FontName=font)
xlabel(xlbl,FontName=font)
ylabel(ylbl,FontName=font)
ylim(limY(2,:))
grid
hold on
for j=1:length(linee)
    line(linee(j)*[1,1],limY(2,:), 'Color','black')
end
text(linee(:)+0.5*[1,1,1,1/2,1]',textLimY(2,:),str,FontSize=6.5, FontWeight="bold",FontName=font)
% text(linee(:)+0.5*[1,1,1,1]',textLimY(2,:),str,FontSize=6.5, FontWeight="bold",FontName=font)

subplot(3,1,3)
plot(t,fun(:,3),LineWidth=1,color="b")
subtitle(axs(3),FontName=font)
xlabel(xlbl,FontName=font)
ylabel(ylbl,FontName=font)
ylim(limY(3,:))
grid
hold on
for j=1:length(linee)
    line(linee(j)*[1,1],limY(3,:), 'Color','black')
end
text(linee(:)+0.5*[1,1,1,1/2,1]',textLimY(3,:),str,FontSize=6.5, FontWeight="bold",FontName=font)
% text(linee(:)+0.5*[1,1,1,1]',textLimY(3,:),str,FontSize=6.5, FontWeight="bold",FontName=font)

exportgraphics(f,"..\Relazione\4_Dati\img\"+"Acc LungaF"+".pdf","ContentType","vector")


