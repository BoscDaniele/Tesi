clear 
close all
clc

%% Rilevamento
path="dbdm\salita-discesa2\";
rilievo=2; % il rilievo 0 è la gravità, il rilievo 1 è per stabilire di quanto è ruotato il sensore rispetto all'asse z


%% Sistema di riferimento sensore = sistema di riferimento bicicletta

[gzRot,gMedio] = GZRot(path);


