close all
clear
clc

inizio=2;

dbA=importdata("dbdm\nuovo\BlueCoin_Log_N001.csv");
acc=dbA.data(inizio:end,2:4);
t=dbA.data(inizio:end,1)*1e-3;
t=t-t(1);

% acc=lowpass(acc,40,1e3);
% acc=highpass(acc,1,1e1);

figure(Name="Data (t)")
subplot(3,1,1)
plot(t,acc(:,1),LineWidth=1);
grid;
title("Data(t)");
subtitle("X");
subplot(3,1,2)
plot(t,acc(:,2),LineWidth=1);
grid;
subtitle("Y");
subplot(3,1,3)
plot(t,acc(:,3),LineWidth=1);
grid;
subtitle("Z");


%% Dft
facc=fft(acc);
f=1/t(2);
periodo=t(2);
l=length(t);

figure
plot(f/l*(0:l-1),abs(facc(:,1)),LineWidth=1);
grid;
title("Data(f)");
subtitle("X");

figure(Name="Data (f)")
subplot(3,1,1)
plot(f/l*(0:l-1),facc(:,1),LineWidth=1);
grid;
title("Data(f)");
subtitle("X");
subplot(3,1,2)
plot(f/l*(0:l-1),abs(facc(:,2)),LineWidth=1);
grid;
subtitle("Y");
subplot(3,1,3)
plot(f/l*(0:l-1),abs(facc(:,3)),LineWidth=1);
grid;
subtitle("Z");