clear all
close all
clc

% lunga
path=".\dati\lunga_forte\";
% utili=[2,4,6];
rilievo=2;

% % curvaU
% path=".\dati\curvaU_forte\";
% % utili: 2,4,6
% rilievo=;

[gzRot,gMedio] = GZRot(path);

sr = 25; %sample rate

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*gzRot*9.81/-gMedio;
vang_rotta=(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1);

[t,acc,vang,mag]=AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta);

filtered_acc=lowpass(acc,0.5,sr);
mean_acc=movmean(acc,40);


% c=["r","g","b"];
% str="Accelerazione";
% sub=["asse X", "asse Y", "asse Z"];
%
% lvl=[0,10,20];
%
% figure("Name",str)
% for i=1:3
%     subplot(3,1,i)
%     if(i==1)
%         title(str)
%     end
%     plot(t,acc(:,i),LineWidth=1,Color=c(i))
%     subtitle(sub(i));
%     xlabel("t(s)")
%     ylabel("m/s^2")
%     grid
%     % hold on
%     % plot(tF,vangF(:,i),LineWidth=1)
%     % plot(tF,ones(length(tF))*lvl(i),LineWidth=1,Color="black")
% end


c=["r","g","b"];
str="Accelerazione";
sub=["asse X", "asse Y", "asse Z"];

lvl=[0,10,20];

figure("Name",str)
subplot(3,1,1)
title(str)
plot(t,mean_acc(:,1),LineWidth=1,Color=c(1))
subtitle(sub(1));
xlabel("t(s)")
ylabel("m/s^2")
% ylim([-10,10])
grid
subplot(3,1,2)
plot(t,mean_acc(:,2),LineWidth=1,Color=c(2))
subtitle(sub(2));
xlabel("t(s)")
ylabel("m/s^2")
% ylim([-5,5])
grid
subplot(3,1,3)
plot(t,mean_acc(:,3),LineWidth=1,Color=c(3))
subtitle(sub(3));
xlabel("t(s)")
ylabel("m/s^2")
grid
% hold on
% plot(tF,vangF(:,i),LineWidth=1)
% plot(tF,ones(length(tF))*lvl(i),LineWidth=1,Color="black")




c=["r","g","b"];
str="Accelerazione Filtrata";
sub=["asse X", "asse Y", "asse Z"];

lvl=[0,10,20];

figure("Name",str)
for i=1:3
    subplot(3,1,i)
    if(i==1)
        title(str)
    end
    plot(t,filtered_acc(:,i),LineWidth=1,Color=c(i))
    subtitle(sub(i));
    xlabel("t(s)")
    ylabel("m/s^2")
    grid
    % hold on
    % plot(tF,vangF(:,i),LineWidth=1)
    % plot(tF,ones(length(tF))*lvl(i),LineWidth=1,Color="black")
end


c=["r","g","b"];
str="Velocità Angolare";
sub=["Roll", "Pitch", "Yaw"];

lvl=[0,10,20];

figure("Name",str)
for i=1:3
    subplot(3,1,i)
    if(i==1)
        title(str)
    end
    plot(t,vang(:,i),LineWidth=1,Color=c(i))
    subtitle(sub(i));
    xlabel("t(s)")
    ylabel("deg/s")
    grid
    % hold on
    % plot(tF,vangF(:,i),LineWidth=1)
    % plot(tF,ones(length(tF))*lvl(i),LineWidth=1,Color="black")
end


c=["r","g","b"];
str="Campo Magnetico";
sub=["asse X", "asse Y", "asse Z"];

lvl=[0,10,20];

figure("Name",str)
for i=1:3
    subplot(3,1,i)
    if(i==1)
        title(str)
    end
    plot(t,mag(:,i),LineWidth=1,Color=c(i))
    subtitle(sub(i));
    xlabel("t(s)")
    ylabel("µT")
    grid
    % hold on
    % plot(tF,vangF(:,i),LineWidth=1)
    % plot(tF,ones(length(tF))*lvl(i),LineWidth=1,Color="black")
end

dim1=[0.15 0.2 0.284 0.6];
dim2=[0.436 0.24 0.09 0.525];

figure
plot(t,acc(:,1),LineWidth=1,Color="r");
xlabel("t(s)")
ylabel("m/s^2")
ylim([-10,10])
grid
hold on
plot(t,mean_acc(:,1),LineWidth=1,Color="black")
legend("Accelerazione","Accelerazione Media")
annotation("rectangle",dim1,"Color","b",LineWidth=1.5)
annotation("rectangle",dim2,"Color","g",LineWidth=1.5)
