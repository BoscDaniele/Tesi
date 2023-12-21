clear
close all
clc

%% Rilevamento
path="dbdm\salita-discesa\";
rilievo=1;


%% Sistema di riferimento sensore = sistema di riferimento bicicletta
inizioG=1;
fineG=inizioG+675;

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;
[gRot,gMedio]=GRot(inizioG,fineG,db);

%% Impostazioni
inizio=1;
fine=length(db);

%% Import e filtraggio dati
t=db(inizio:fine,1)*1e-3;
t=t-t(1);

acc=db(inizio:fine,2:4)*gRot*9.81/-gMedio;
vang=db(inizio:fine,5:7)*gRot*2*pi/360*1e-3;

modelname = "provaSimulink";
open_system(modelname)

set_param(modelname,"StartTime","0.04","StopTime",num2str(t(end)))

sim(modelname);
