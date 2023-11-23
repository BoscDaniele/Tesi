close all
clear
clc

%import dati accelerazioni
dbA=importdata("dbdm\Acceleration.csv");
acc=dbA.data;

%estrazione tempi
t=str2num(cell2mat(dbA.textdata(4:end,4)));
t=t-t(1);

%plot
figure(Name="Accelereazione")
subplot(3,1,1);
plot(t,acc(:,1),LineWidth=1,Color="b");
subtitle("X", Color="b");
subplot(3,1,2);
plot(t,acc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="b");
subplot(3,1,3);
plot(t,acc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="b");

%import dati accelerazioni
dbA1=importdata("dbdm\Acceleration (1).csv");
acc1=dbA1.data;

%estrazione tempi
t1=str2num(cell2mat(dbA1.textdata(4:end,4)));
t1=t1-t1(1);

%plot
figure(Name="Accelereazione 1")
subplot(3,1,1);
plot(t1,acc1(:,1),LineWidth=1,Color="b");
subtitle("X", Color="b");
subplot(3,1,2);
plot(t1,acc1(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="b");
subplot(3,1,3);
plot(t1,acc1(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="b");

%import dati accelerazioni
dbA2=importdata("dbdm\Acceleration (2).csv");
acc2=dbA2.data;

%estrazione tempi
t2=str2num(cell2mat(dbA2.textdata(4:end,4)));
t2=t2-t2(1);

%plot
figure(Name="Accelereazione 2")
subplot(3,1,1);
plot(t2,acc2(:,1),LineWidth=1,Color="b");
subtitle("X", Color="b");
subplot(3,1,2);
plot(t2,acc2(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="b");
subplot(3,1,3);
plot(t2,acc2(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="b");
