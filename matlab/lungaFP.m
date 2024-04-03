clear
close all
clc

%% Import dati

% lunga
path=".\dati\lunga";

% i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% rotazione per far coincidere il sistema di riferimento dell'accelerometro
% con quello della bicicletta
rilievoP=5;
rilievoF=4;

% % curva
% path=".\dati\curva";
% 
% % i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% % rotazione per far coincidere il sistema di riferimento dell'accelerometro
% % con quello della bicicletta
% rilievoP=4;
% rilievoF=2;

% % curvaU
% path=".\dati\curvaU";
% 
% % i rilievi 0 e 1 vengono utilizzati solo per il calcolo della matrice di
% % rotazione per far coincidere il sistema di riferimento dell'accelerometro
% % con quello della bicicletta
% rilievoP=3;
% rilievoF=4;

[gzRotP,gMedioP] = GZRot(path+"_piano\");
[gzRotF,gMedioF] = GZRot(path+"_forte\");

sr = 25; %sample rate


%% Import Dati
dbP=importdata(path + "_piano\" + "BlueCoin_Log_N00"+rilievoP+".csv").data;
dbF=importdata(path + "_forte\" + "BlueCoin_Log_N00"+rilievoF+".csv").data;

%% Estrazione e Correzione Dati Piano
% estrazione dati tempo e conversione in secondi
tP_rotta=dbP(:,1)*1e-3;
tP_rotta=tP_rotta-tP_rotta(1);

accP_rotta=dbP(:,2:4)*gzRotP*9.81/-gMedioP;
vangP_rotta=(dbP(:,5:7)*1e-3);
magP_rotta=([dbP(:,8),-dbP(:,9),dbP(:,10)]*1e-1);


tP=[tP_rotta(1)];
lastP=1;

accP=[accP_rotta(1,:)];
vangP=[vangP_rotta(1,:)];
magP=[magP_rotta(1,:)];

for i=2:length(tP_rotta)
    % disp(num2str(tP(i)-tempoP(end)))
    if(tP_rotta(i)-tP(end)>=0.039 && tP_rotta(i)-tP(end)<0.041)
        tP=[tP,tP_rotta(i)];
        accP=[accP;accP_rotta(i,:)];
        vangP=[vangP;vangP_rotta(i,:)];
        magP=[magP;magP_rotta(i,:)];
        lastP=i;

        if((i-lastP)>5)
            disp("lammerda, abbiamo saltato 0.004s "+num2str(i));
            brake;
        end
    end
   
end


diff_tP=zeros(length(tP),1);
diff_tP(1)=tP(1);
for i=2:length(tP)
    diff_tP(i)=tP(i)-tP(i-1);
end


disp("Piano: max diff_t: "+num2str(max(diff_tP(2:end))));
disp("Piano: min diff_t: "+num2str(min(diff_tP(2:end))));


%% Estrazione e Correzione Dati Forte
tF_rotta=dbF(:,1)*1e-3;
tF_rotta=tF_rotta-tF_rotta(1);

accF_rotta=dbF(:,2:4)*gzRotF*9.81/-gMedioF;
vangF_rotta=(dbF(:,5:7)*1e-3);
magF_rotta=([dbF(:,8),-dbF(:,9),dbF(:,10)]*1e-1);

tF=[tF_rotta(1)];
lastF=1;

accF=[accF_rotta(1,:)];
vangF=[vangF_rotta(1,:)];
magF=[magF_rotta(1,:)];

for i=2:length(tF_rotta)
    % disp(num2str(tP(i)-tempoP(end)))
    if(tF_rotta(i)-tF(end)>=0.039 && tF_rotta(i)-tF(end)<0.041)
        tF=[tF,tF_rotta(i)];
        accF=[accF;accF_rotta(i,:)];
        vangF=[vangF;vangF_rotta(i,:)];
        magF=[magF;magF_rotta(i,:)];
        lastF=i;

        if((i-lastF)>5)
            disp("lammerda, abbiamo saltato 0.004s "+num2str(i));
            brake;
        end
    end
   
end


diff_tF=zeros(length(tF),1);
diff_tF(1)=tF(1);
for i=2:length(tF)
    diff_tF(i)=tF(i)-tF(i-1);
end

disp("Forte: max diff_t: "+num2str(max(diff_tF(2:end))));
disp("Forte: min diff_t: "+num2str(min(diff_tF(2:end))));


%% Figure

% Accelerazione
acc_fig=1;
acc_mean_fig=1; %1
acc_median_fig=0;

acc_var_fig=1; %1
acc_std_fig=1; %1
acc_rms_fig=1; %1
acc_kurt_fig=0;

acc_norm_fig=0;
accXY_norm_fig=0;

acc_max_fig=1; %1
acc_min_fig=0;
acc_peak_fig=1; %1

% Accelerazione Filtrata
filtered_acc_fig=1;
filtered_acc_mean_fig=0;
filtered_acc_median_fig=0;

filtered_acc_var_fig=0; %1
filtered_acc_std_fig=0; %1
filtered_acc_rms_fig=0; %1
filtered_acc_kurt_fig=0;

filtered_acc_norm_fig=0; %1
filtered_accXY_norm_fig=0; %1

filtered_acc_max_fig=0; %1
filtered_acc_min_fig=0;
filtered_acc_peak_fig=0; %1

% Velocità
vel_fig=1;
vel_mean_fig=0;
vel_median_fig=0;

vel_var_fig=0;
vel_std_fig=0;
vel_rms_fig=0;
vel_kurt_fig=0;

vel_norm_fig=0;
velXY_norm_fig=0; %1

vel_max_fig=0;
vel_min_fig=0;
vel_peak_fig=0; %1

% Velocità Angolare
vang_fig=0;
vang_mean_fig=0; %1
vang_median_fig=0;

vang_var_fig=0; %1
vang_std_fig=0; %1
vang_rms_fig=0; %1
vang_kurt_fig=0;

vang_norm_fig=0;
vangXY_norm_fig=0;

vang_max_fig=0; %1
vang_min_fig=0; %1
vang_peak_fig=0;

% Angoli
ang_fig=0;
ang_mean_fig=0;
ang_median_fig=0;

ang_var_fig=0;
ang_std_fig=0;
ang_rms_fig=0;
ang_kurt_fig=0;

ang_norm_fig=0;
angXY_norm_fig=0;

ang_max_fig=0;
ang_min_fig=0;
ang_peak_fig=0;

% Campo Magnetico
mag_fig=0;
mag_mean_fig=0;
mag_median_fig=0;

mag_var_fig=0;
mag_std_fig=0;
mag_rms_fig=0;
mag_kurt_fig=0;

mag_norm_fig=0;
magXY_norm_fig=0;

mag_max_fig=0;
mag_min_fig=0;
mag_peak_fig=0;


%% Accelerazione
if(acc_fig)

    % Accelerazione Media
    if(acc_mean_fig)
        n_mean=[40,0];
        accP_mean=movmean(accP,n_mean);
        accF_mean=movmean(accF,n_mean);
    end

    % Accelerazione Mediana
    if(acc_median_fig)
        n_median=[40,0];
        accP_median=movmedian(accP,n_median);
        accF_median=movmedian(accF,n_median);
    end

    % Accelerazione Varianza
    if(acc_var_fig)
        n_var=[40,0];
        accP_var=movvar(accP,n_var);
        accF_var=movvar(accF,n_var);
    end

    % Accelerazione Deviazione Standard
    if(acc_std_fig)
        n_std=[40,0];
        accP_std=movstd(accP,n_std);
        accF_std=movstd(accF,n_std);
    end

    % Accelerazione Scarto Quadratico Medio
    if(acc_rms_fig)
        n_rms=[40,0];
        accP_rms=movrms(accP,n_rms);
        accF_rms=movrms(accF,n_rms);
    end

    % Accelerazione Kurtosi
    if(acc_kurt_fig)
        n_kurt=[40,0];
        accP_kurt=movkurt(accP,n_kurt);
        accF_kurt=movkurt(accF,n_kurt);
    end

    % Accelerazione Norma
    if(acc_norm_fig)
        accP_norm=zeros(length(accP),1);
        accF_norm=zeros(length(accF),1);

        for i=1:length(accP)
            accP_norm(i)=norm(accP(i,:));
        end

        for j=1:length(accF)
            accF_norm(j)=norm(accF(j,:));
        end

    end

    % Accelerazione XY Norma
    if(accXY_norm_fig)
        accP_XY_norm=zeros(length(accP),1);
        accF_XY_norm=zeros(length(accF),1);

        for i=1:length(accP)
            accP_XY_norm(i)=norm(accP(i,1:2));
        end

        for j=1:length(accF)
            accF_XY_norm(j)=norm(accF(j,1:2));
        end

    end

    % Accelerazione Max
    if(acc_max_fig || acc_peak_fig)
        n_max=[10,0];
        accP_max=movmax(accP,n_max);
        accF_max=movmax(accF,n_max);
    end

    % Accelerazione Min
    if(acc_min_fig || acc_peak_fig)
        n_min=[10,0];
        accP_min=movmin(accP,n_min);
        accF_min=movmin(accF,n_min);
    end

    % Accelerazione Peak2Peak
    if(acc_peak_fig)
        accP_peak=accP_max-accP_min;
        accF_peak=accF_max-accF_min;
    end

end

%% Accelerazione Filtrata
filtered_accP=lowpass(accP,0.5,sr);
filtered_accF=lowpass(accF,0.5,sr);

if(filtered_acc_fig)

    % Accelerazione Filtrata Media
    if(filtered_acc_mean_fig)
        n_mean=[40,0];
        filtered_accP_mean=movmean(filtered_accP,n_mean);
        filtered_accF_mean=movmean(filtered_accF,n_mean);
    end

    % Accelerazione Filtrata Mediana
    if(filtered_acc_median_fig)
        n_median=[40,0];
        filtered_accP_median=movmedian(filtered_accP,n_median);
        filtered_accF_median=movmedian(filtered_accF,n_median);
    end

    % Accelerazione Filtrata Varianza
    if(filtered_acc_var_fig)
        n_var=[40,0];
        filtered_accP_var=movvar(filtered_accP,n_var);
        filtered_accF_var=movvar(filtered_accF,n_var);
    end

    % Accelerazione Filtrata Deviazione Standard
    if(filtered_acc_std_fig)
        n_std=[40,0];
        filtered_accP_std=movstd(filtered_accP,n_std);
        filtered_accF_std=movstd(filtered_accF,n_std);
    end

    % Accelerazione Filtrata Scarto Quadratico Medio
    if(filtered_acc_rms_fig)
        n_rms=[40,0];
        filtered_accP_rms=movrms(filtered_accP,n_rms);
        filtered_accF_rms=movrms(filtered_accF,n_rms);
    end

    % Accelerazione Filtrata Kurtosi
    if(filtered_acc_kurt_fig)
        n_kurt=[40,0];
        filtered_accP_kurt=movkurt(filtered_accP,n_kurt);
        filtered_accF_kurt=movkurt(filtered_accF,n_kurt);
    end

    % Accelerazione Filtrata Norma
    if(filtered_acc_norm_fig)
        filtered_accP_norm=zeros(length(filtered_accP),1);
        filtered_accF_norm=zeros(length(filtered_accF),1);

        for i=1:length(filtered_accP)
            filtered_accP_norm(i)=norm(filtered_accP(i,:));
        end

        for j=1:length(filtered_accF)
            filtered_accF_norm(j)=norm(filtered_accF(j,:));
        end

    end

    % Accelerazione Filtrata XY Norma
    if(filtered_accXY_norm_fig)
        filtered_accP_XY_norm=zeros(length(filtered_accP),1);
        filtered_accF_XY_norm=zeros(length(filtered_accF),1);

        for i=1:length(filtered_accP)
            filtered_accP_XY_norm(i)=norm(filtered_accP(i,1:2));
        end

        for j=1:length(filtered_accF)
            filtered_accF_XY_norm(j)=norm(filtered_accF(j,1:2));
        end

    end

    % Accelerazione Filtrata Max
    if(filtered_acc_max_fig || filtered_acc_peak_fig)
        n_max=[10,0];
        filtered_accP_max=movmax(filtered_accP,n_max);
        filtered_accF_max=movmax(filtered_accF,n_max);
    end

    % Accelerazione Filtrata Min
    if(filtered_acc_min_fig || filtered_acc_peak_fig)
        n_min=[10,0];
        filtered_accP_min=movmin(filtered_accP,n_min);
        filtered_accF_min=movmin(filtered_accF,n_min);
    end

    % Accelerazione Filtrata Peak2Peak
    if(filtered_acc_peak_fig)
        filtered_accP_peak=filtered_accP_max-filtered_accP_min;
        filtered_accF_peak=filtered_accF_max-filtered_accF_min;
    end

end

%% Velocità
velP=cumsum(accP)*0.04;
velF=cumsum(accF)*0.04;

if(vel_fig)

    % Velocità Media
    if(vel_mean_fig)
        n_mean=[40,0];
        velP_mean=movmean(velP,n_mean);
        velF_mean=movmean(velF,n_mean);
    end

    % Velocità Mediana
    if(vel_median_fig)
        n_median=[40,0];
        velP_median=movmedian(velP,n_median);
        velF_median=movmedian(velF,n_median);
    end

    % Velocità Varianza
    if(vel_var_fig)
        n_var=[40,0];
        velP_var=movvar(velP,n_var);
        velF_var=movvar(velF,n_var);
    end

    % Velocità Deviazione Standard
    if(vel_std_fig)
        n_std=[40,0];
        velP_std=movstd(velP,n_std);
        velF_std=movstd(velF,n_std);
    end

    % Velocità Scarto Quadratico Medio
    if(vel_rms_fig)
        n_rms=[40,0];
        velP_rms=movrms(velP,n_rms);
        velF_rms=movrms(velF,n_rms);
    end

    % Velocità Kurtosi
    if(vel_kurt_fig)
        n_kurt=[40,0];
        velP_kurt=movkurt(velP,n_kurt);
        velF_kurt=movkurt(velF,n_kurt);
    end

    % Velocità Norma
    if(vel_norm_fig)
        velP_norm=zeros(length(velP),1);
        velF_norm=zeros(length(velF),1);

        for i=1:length(velP)
            velP_norm(i)=norm(velP(i,:));
        end

        for j=1:length(velF)
            velF_norm(j)=norm(velF(j,:));
        end

    end

    % Velocità XY Norma
    if(velXY_norm_fig)
        velP_XY_norm=zeros(length(velP),1);
        velF_XY_norm=zeros(length(velF),1);

        for i=1:length(velP)
            velP_XY_norm(i)=norm(velP(i,1:2));
        end

        for j=1:length(velF)
            velF_XY_norm(j)=norm(velF(j,1:2));
        end

    end

    % Velocità Max
    if(vel_max_fig || vel_peak_fig)
        n_max=[10,0];
        velP_max=movmax(velP,n_max);
        velF_max=movmax(velF,n_max);
    end

    % Velocità Min
    if(vel_min_fig || vel_peak_fig)
        n_min=[10,0];
        velP_min=movmin(velP,n_min);
        velF_min=movmin(velF,n_min);
    end

    % Velocità Peak2Peak
    if(vel_peak_fig)
        velP_peak=velP_max-velP_min;
        velF_peak=velF_max-velF_min;
    end

end


%% Velocità Angolare

if(vang_fig)

    % Velocità Angolare Media
    if(vang_mean_fig)
        n_mean=[40,0];
        vangP_mean=movmean(vangP,n_mean);
        vangF_mean=movmean(vangF,n_mean);
    end

    % Velocità Angolare Mediana
    if(vang_median_fig)
        n_median=[40,0];
        vangP_median=movmedian(vangP,n_median);
        vangF_median=movmedian(vangF,n_median);
    end

    % Velocità Angolare Varianza
    if(vang_var_fig)
        n_var=[40,0];
        vangP_var=movvar(vangP,n_var);
        vangF_var=movvar(vangF,n_var);
    end

    % Velocità Angolare Deviazione Standard
    if(vang_std_fig)
        n_std=[40,0];
        vangP_std=movstd(vangP,n_std);
        vangF_std=movstd(vangF,n_std);
    end

    % Velocità Angolare Scarto Quadratico Medio
    if(vang_rms_fig)
        n_rms=[40,0];
        vangP_rms=movrms(vangP,n_rms);
        vangF_rms=movrms(vangF,n_rms);
    end

    % Velocità Angolare Kurtosi
    if(vang_kurt_fig)
        n_kurt=[40,0];
        vangP_kurt=movkurt(vangP,n_kurt);
        vangF_kurt=movkurt(vangF,n_kurt);
    end

    % Velocità Angolare Norma
    if(vang_norm_fig)
        vangP_norm=zeros(length(vangP),1);
        vangF_norm=zeros(length(vangF),1);

        for i=1:length(vangP)
            vangP_norm(i)=norm(vangP(i,:));
        end

        for j=1:length(vangF)
            vangF_norm(j)=norm(vangF(j,:));
        end

    end

    % Velocità Angolare XY Norma
    if(vangXY_norm_fig)
        vangP_XY_norm=zeros(length(vangP),1);
        vangF_XY_norm=zeros(length(vangF),1);

        for i=1:length(vangP)
            vangP_XY_norm(i)=norm(vangP(i,1:2));
        end

        for j=1:length(vangF)
            vangF_XY_norm(j)=norm(vangF(j,1:2));
        end

    end

    % Velocità Angolare Max
    if(vang_max_fig || vang_peak_fig)
        n_max=[10,0];
        vangP_max=movmax(vangP,n_max);
        vangF_max=movmax(vangF,n_max);
    end

    % Velocità Angolare Min
    if(vang_min_fig || vang_peak_fig)
        n_min=[10,0];
        vangP_min=movmin(vangP,n_min);
        vangF_min=movmin(vangF,n_min);
    end

    % Velocità Angolare Peak2Peak
    if(vang_peak_fig)
        vangP_peak=vangP_max-vangP_min;
        vangF_peak=vangF_max-vangF_min;
    end

end

%% Angoli
angP=cumsum(vangP)*0.04;
angF=cumsum(vangF)*0.04;

if(ang_fig)

    % Angoli Media
    if(ang_mean_fig)
        n_mean=[40,0];
        angP_mean=movmean(angP,n_mean);
        angF_mean=movmean(angF,n_mean);
    end

    % Angoli Mediana
    if(ang_median_fig)
        n_median=[40,0];
        angP_median=movmedian(angP,n_median);
        angF_median=movmedian(angF,n_median);
    end

    % Angoli Varianza
    if(ang_var_fig)
        n_var=[40,0];
        angP_var=movvar(angP,n_var);
        angF_var=movvar(angF,n_var);
    end

    % Angoli Deviazione Standard
    if(ang_std_fig)
        n_std=[40,0];
        angP_std=movstd(angP,n_std);
        angF_std=movstd(angF,n_std);
    end

    % Angoli Scarto Quadratico Medio
    if(ang_rms_fig)
        n_rms=[40,0];
        angP_rms=movrms(angP,n_rms);
        angF_rms=movrms(angF,n_rms);
    end

    % Angoli Kurtosi
    if(ang_kurt_fig)
        n_kurt=[40,0];
        angP_kurt=movkurt(angP,n_kurt);
        angF_kurt=movkurt(angF,n_kurt);
    end

    % Angoli Norma
    if(ang_norm_fig)
        angP_norm=zeros(length(angP),1);
        angF_norm=zeros(length(angF),1);

        for i=1:length(angP)
            angP_norm(i)=norm(angP(i,:));
        end

        for j=1:length(angF)
            angF_norm(j)=norm(angF(j,:));
        end

    end

    % Angoli XY Norma
    if(angXY_norm_fig)
        angP_XY_norm=zeros(length(angP),1);
        angF_XY_norm=zeros(length(angF),1);

        for i=1:length(angP)
            angP_XY_norm(i)=norm(angP(i,1:2));
        end

        for j=1:length(angF)
            angF_XY_norm(j)=norm(angF(j,1:2));
        end

    end

    % Angoli Max
    if(ang_max_fig || ang_peak_fig)
        n_max=[10,0];
        angP_max=movmax(angP,n_max);
        angF_max=movmax(angF,n_max);
    end

    % Angoli Min
    if(ang_min_fig || ang_peak_fig)
        n_min=[10,0];
        angP_min=movmin(angP,n_min);
        angF_min=movmin(angF,n_min);
    end

    % Angoli Peak2Peak
    if(ang_peak_fig)
        angP_peak=angP_max-angP_min;
        angF_peak=angF_max-angF_min;
    end

end

%% Campo Magnetico

if(mag_fig)

    % Campo Magnetico Media
    if(mag_mean_fig)
        n_mean=[40,0];
        magP_mean=movmean(magP,n_mean);
        magF_mean=movmean(magF,n_mean);
    end

    % Campo Magnetico Mediana
    if(mag_median_fig)
        n_median=[40,0];
        magP_median=movmedian(magP,n_median);
        magF_median=movmedian(magF,n_median);
    end

    % Campo Magnetico Varianza
    if(mag_var_fig)
        n_var=[40,0];
        magP_var=movvar(magP,n_var);
        magF_var=movvar(magF,n_var);
    end

    % Campo Magnetico Deviazione Standard
    if(mag_std_fig)
        n_std=[40,0];
        magP_std=movstd(magP,n_std);
        magF_std=movstd(magF,n_std);
    end

    % Campo Magnetico Scarto Quadratico Medio
    if(mag_rms_fig)
        n_rms=[40,0];
        magP_rms=movrms(magP,n_rms);
        magF_rms=movrms(magF,n_rms);
    end

    % Campo Magnetico Kurtosi
    if(mag_kurt_fig)
        n_kurt=[40,0];
        magP_kurt=movkurt(magP,n_kurt);
        magF_kurt=movkurt(magF,n_kurt);
    end

    % Campo Magnetico Norma
    if(mag_norm_fig)
        magP_norm=zeros(length(magP),1);
        magF_norm=zeros(length(magF),1);

        for i=1:length(magP)
            magP_norm(i)=norm(magP(i,:));
        end

        for j=1:length(magF)
            magF_norm(j)=norm(magF(j,:));
        end

    end

    % Campo Magnetico XY Norma
    if(magXY_norm_fig)
        magP_XY_norm=zeros(length(magP),1);
        magF_XY_norm=zeros(length(magF),1);

        for i=1:length(magP)
            magP_XY_norm(i)=norm(magP(i,1:2));
        end

        for j=1:length(magF)
            magF_XY_norm(j)=norm(magF(j,1:2));
        end

    end

    % Campo Magnetico Max
    if(mag_max_fig || mag_peak_fig)
        n_max=[10,0];
        magP_max=movmax(magP,n_max);
        magF_max=movmax(magF,n_max);
    end

    % Campo Magnetico Min
    if(mag_min_fig || mag_peak_fig)
        n_min=[10,0];
        magP_min=movmin(magP,n_min);
        magF_min=movmin(magF,n_min);
    end

    % Campo Magnetico Peak2Peak
    if(mag_peak_fig)
        magP_peak=magP_max-magP_min;
        magF_peak=magF_max-magF_min;
    end

end



%% Figure
lower_lim=zeros(3,1);
upper_lim=zeros(3,1);

% Accelerazione
str="Accelerazione";
unit="m/s^2";

if(acc_fig)
    for i=1:3
        if(min(accP(:,i))<min(accF(:,i)))
            lower_lim(i)=floor(min(accP(:,i)));
        else
            lower_lim(i)=floor(min(accF(:,i)));
        end

        if(max(accP(:,i))>max(accF(:,i)))
            upper_lim(i)=ceil(max(accP(:,i)));
        else
            upper_lim(i)=ceil(max(accF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,accP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,accP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,accP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,accF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,accF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,accF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid


    if(acc_mean_fig)
        for i=1:3
            if(min(accP_mean(:,i))<min(accF_mean(:,i)))
                lower_lim(i)=floor(min(accP_mean(:,i)));
            else
                lower_lim(i)=floor(min(accF_mean(:,i)));
            end

            if(max(accP_mean(:,i))>max(accF_mean(:,i)))
                upper_lim(i)=ceil(max(accP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(accF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,accP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,accP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_median_fig)
        for i=1:3
            if(min(accP_median(:,i))<min(accF_median(:,i)))
                lower_lim(i)=floor(min(accP_median(:,i)));
            else
                lower_lim(i)=floor(min(accF_median(:,i)));
            end

            if(max(accP_median(:,i))>max(accF_median(:,i)))
                upper_lim(i)=ceil(max(accP_median(:,i)));
            else
                upper_lim(i)=ceil(max(accF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,accP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,accP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(acc_var_fig)
        for i=1:3
            if(min(accP_var(:,i))<min(accF_var(:,i)))
                lower_lim(i)=floor(min(accP_var(:,i)));
            else
                lower_lim(i)=floor(min(accF_var(:,i)));
            end

            if(max(accP_var(:,i))>max(accF_var(:,i)))
                upper_lim(i)=ceil(max(accP_var(:,i)));
            else
                upper_lim(i)=ceil(max(accF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,accP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,accP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_std_fig)
        for i=1:3
            if(min(accP_std(:,i))<min(accF_std(:,i)))
                lower_lim(i)=floor(min(accP_std(:,i)));
            else
                lower_lim(i)=floor(min(accF_std(:,i)));
            end

            if(max(accP_std(:,i))>max(accF_std(:,i)))
                upper_lim(i)=ceil(max(accP_std(:,i)));
            else
                upper_lim(i)=ceil(max(accF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,accP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,accP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_rms_fig)
        for i=1:3
            if(min(accP_rms(:,i))<min(accF_rms(:,i)))
                lower_lim(i)=floor(min(accP_rms(:,i)));
            else
                lower_lim(i)=floor(min(accF_rms(:,i)));
            end

            if(max(accP_rms(:,i))>max(accF_rms(:,i)))
                upper_lim(i)=ceil(max(accP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(accF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,accP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,accP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_kurt_fig)
        for i=1:3
            if(min(accP_kurt(:,i))<min(accF_kurt(:,i)))
                lower_lim(i)=floor(min(accP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(accF_kurt(:,i)));
            end

            if(max(accP_kurt(:,i))>max(accF_kurt(:,i)))
                upper_lim(i)=ceil(max(accP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(accF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,accP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,accP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_norm_fig)
        if(min(accP_norm)<min(accF_norm))
            lower_lim(1)=floor(min(accP_norm));
        else
            lower_lim(1)=floor(min(accF_norm));
        end

        if(max(accP_norm)>max(accF_norm))
            upper_lim(1)=ceil(max(accP_norm));
        else
            upper_lim(1)=ceil(max(accF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,accP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,accP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,accF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(accXY_norm_fig)
        if(min(accP_XY_norm)<min(accF_XY_norm))
            lower_lim(1)=floor(min(accP_XY_norm));
        else
            lower_lim(1)=floor(min(accF_XY_norm));
        end

        if(max(accP_XY_norm)>max(accF_XY_norm))
            upper_lim(1)=ceil(max(accP_XY_norm));
        else
            upper_lim(1)=ceil(max(accF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,accP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,accP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,accF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(acc_max_fig)
        for i=1:3
            if(min(accP_max(:,i))<min(accF_max(:,i)))
                lower_lim(i)=floor(min(accP_max(:,i)));
            else
                lower_lim(i)=floor(min(accF_max(:,i)));
            end

            if(max(accP_max(:,i))>max(accF_max(:,i)))
                upper_lim(i)=ceil(max(accP_max(:,i)));
            else
                upper_lim(i)=ceil(max(accF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,accP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,accP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_min_fig)
        for i=1:3
            if(min(accP_min(:,i))<min(accF_min(:,i)))
                lower_lim(i)=floor(min(accP_min(:,i)));
            else
                lower_lim(i)=floor(min(accF_min(:,i)));
            end

            if(max(accP_min(:,i))>max(accF_min(:,i)))
                upper_lim(i)=ceil(max(accP_min(:,i)));
            else
                upper_lim(i)=ceil(max(accF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,accP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,accP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(acc_peak_fig)
        for i=1:3
            if(min(accP_peak(:,i))<min(accF_peak(:,i)))
                lower_lim(i)=floor(min(accP_peak(:,i)));
            else
                lower_lim(i)=floor(min(accF_peak(:,i)));
            end

            if(max(accP_peak(:,i))>max(accF_peak(:,i)))
                upper_lim(i)=ceil(max(accP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(accF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,accP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,accP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,accP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,accF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,accP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,accP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,accP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,accF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,accF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,accF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end

% Accelerazione Filtrata
str="Accelerazione Filtrata";
unit="m/s^2";

if(filtered_acc_fig)
    for i=1:3
        if(min(filtered_accP(:,i))<min(filtered_accF(:,i)))
            lower_lim(i)=floor(min(filtered_accP(:,i)));
        else
            lower_lim(i)=floor(min(filtered_accF(:,i)));
        end

        if(max(filtered_accP(:,i))>max(filtered_accF(:,i)))
            upper_lim(i)=ceil(max(filtered_accP(:,i)));
        else
            upper_lim(i)=ceil(max(filtered_accF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,filtered_accP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,filtered_accP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,filtered_accP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,filtered_accF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,filtered_accF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,filtered_accF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid


    if(filtered_acc_mean_fig)
        for i=1:3
            if(min(filtered_accP_mean(:,i))<min(filtered_accF_mean(:,i)))
                lower_lim(i)=floor(min(filtered_accP_mean(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_mean(:,i)));
            end

            if(max(filtered_accP_mean(:,i))>max(filtered_accF_mean(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,filtered_accP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,filtered_accP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_median_fig)
        for i=1:3
            if(min(filtered_accP_median(:,i))<min(filtered_accF_median(:,i)))
                lower_lim(i)=floor(min(filtered_accP_median(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_median(:,i)));
            end

            if(max(filtered_accP_median(:,i))>max(filtered_accF_median(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_median(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,filtered_accP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,filtered_accP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(filtered_acc_var_fig)
        for i=1:3
            if(min(filtered_accP_var(:,i))<min(filtered_accF_var(:,i)))
                lower_lim(i)=floor(min(filtered_accP_var(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_var(:,i)));
            end

            if(max(filtered_accP_var(:,i))>max(filtered_accF_var(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_var(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,filtered_accP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,filtered_accP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_std_fig)
        for i=1:3
            if(min(filtered_accP_std(:,i))<min(filtered_accF_std(:,i)))
                lower_lim(i)=floor(min(filtered_accP_std(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_std(:,i)));
            end

            if(max(filtered_accP_std(:,i))>max(filtered_accF_std(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_std(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,filtered_accP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,filtered_accP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_rms_fig)
        for i=1:3
            if(min(filtered_accP_rms(:,i))<min(filtered_accF_rms(:,i)))
                lower_lim(i)=floor(min(filtered_accP_rms(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_rms(:,i)));
            end

            if(max(filtered_accP_rms(:,i))>max(filtered_accF_rms(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,filtered_accP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,filtered_accP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_kurt_fig)
        for i=1:3
            if(min(filtered_accP_kurt(:,i))<min(filtered_accF_kurt(:,i)))
                lower_lim(i)=floor(min(filtered_accP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_kurt(:,i)));
            end

            if(max(filtered_accP_kurt(:,i))>max(filtered_accF_kurt(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,filtered_accP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,filtered_accP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_norm_fig)
        if(min(filtered_accP_norm)<min(filtered_accF_norm))
            lower_lim(1)=floor(min(filtered_accP_norm));
        else
            lower_lim(1)=floor(min(filtered_accF_norm));
        end

        if(max(filtered_accP_norm)>max(filtered_accF_norm))
            upper_lim(1)=ceil(max(filtered_accP_norm));
        else
            upper_lim(1)=ceil(max(filtered_accF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,filtered_accP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,filtered_accP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,filtered_accF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(filtered_accXY_norm_fig)
        if(min(filtered_accP_XY_norm)<min(filtered_accF_XY_norm))
            lower_lim(1)=floor(min(filtered_accP_XY_norm));
        else
            lower_lim(1)=floor(min(filtered_accF_XY_norm));
        end

        if(max(filtered_accP_XY_norm)>max(filtered_accF_XY_norm))
            upper_lim(1)=ceil(max(filtered_accP_XY_norm));
        else
            upper_lim(1)=ceil(max(filtered_accF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,filtered_accP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,filtered_accP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,filtered_accF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(filtered_acc_max_fig)
        for i=1:3
            if(min(filtered_accP_max(:,i))<min(filtered_accF_max(:,i)))
                lower_lim(i)=floor(min(filtered_accP_max(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_max(:,i)));
            end

            if(max(filtered_accP_max(:,i))>max(filtered_accF_max(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_max(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,filtered_accP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,filtered_accP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_min_fig)
        for i=1:3
            if(min(filtered_accP_min(:,i))<min(filtered_accF_min(:,i)))
                lower_lim(i)=floor(min(filtered_accP_min(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_min(:,i)));
            end

            if(max(filtered_accP_min(:,i))>max(filtered_accF_min(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_min(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,filtered_accP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,filtered_accP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(filtered_acc_peak_fig)
        for i=1:3
            if(min(filtered_accP_peak(:,i))<min(filtered_accF_peak(:,i)))
                lower_lim(i)=floor(min(filtered_accP_peak(:,i)));
            else
                lower_lim(i)=floor(min(filtered_accF_peak(:,i)));
            end

            if(max(filtered_accP_peak(:,i))>max(filtered_accF_peak(:,i)))
                upper_lim(i)=ceil(max(filtered_accP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(filtered_accF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,filtered_accP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,filtered_accP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,filtered_accP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,filtered_accF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,filtered_accP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,filtered_accP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,filtered_accP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,filtered_accF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,filtered_accF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,filtered_accF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end


% Velocità
str="Velocità";
unit="m/s";

if(vel_fig)
    for i=1:3
        if(min(velP(:,i))<min(velF(:,i)))
            lower_lim(i)=floor(min(velP(:,i)));
        else
            lower_lim(i)=floor(min(velF(:,i)));
        end

        if(max(velP(:,i))>max(velF(:,i)))
            upper_lim(i)=ceil(max(velP(:,i)));
        else
            upper_lim(i)=ceil(max(velF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,velP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,velP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,velP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,velF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,velF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,velF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    
    if(vel_mean_fig)
        for i=1:3
            if(min(velP_mean(:,i))<min(velF_mean(:,i)))
                lower_lim(i)=floor(min(velP_mean(:,i)));
            else
                lower_lim(i)=floor(min(velF_mean(:,i)));
            end

            if(max(velP_mean(:,i))>max(velF_mean(:,i)))
                upper_lim(i)=ceil(max(velP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(velF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,velP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,velP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_median_fig)
        for i=1:3
            if(min(velP_median(:,i))<min(velF_median(:,i)))
                lower_lim(i)=floor(min(velP_median(:,i)));
            else
                lower_lim(i)=floor(min(velF_median(:,i)));
            end

            if(max(velP_median(:,i))>max(velF_median(:,i)))
                upper_lim(i)=ceil(max(velP_median(:,i)));
            else
                upper_lim(i)=ceil(max(velF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,velP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,velP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(vel_var_fig)
        for i=1:3
            if(min(velP_var(:,i))<min(velF_var(:,i)))
                lower_lim(i)=floor(min(velP_var(:,i)));
            else
                lower_lim(i)=floor(min(velF_var(:,i)));
            end

            if(max(velP_var(:,i))>max(velF_var(:,i)))
                upper_lim(i)=ceil(max(velP_var(:,i)));
            else
                upper_lim(i)=ceil(max(velF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,velP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,velP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_std_fig)
        for i=1:3
            if(min(velP_std(:,i))<min(velF_std(:,i)))
                lower_lim(i)=floor(min(velP_std(:,i)));
            else
                lower_lim(i)=floor(min(velF_std(:,i)));
            end

            if(max(velP_std(:,i))>max(velF_std(:,i)))
                upper_lim(i)=ceil(max(velP_std(:,i)));
            else
                upper_lim(i)=ceil(max(velF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,velP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,velP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_rms_fig)
        for i=1:3
            if(min(velP_rms(:,i))<min(velF_rms(:,i)))
                lower_lim(i)=floor(min(velP_rms(:,i)));
            else
                lower_lim(i)=floor(min(velF_rms(:,i)));
            end

            if(max(velP_rms(:,i))>max(velF_rms(:,i)))
                upper_lim(i)=ceil(max(velP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(velF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,velP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,velP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_kurt_fig)
        for i=1:3
            if(min(velP_kurt(:,i))<min(velF_kurt(:,i)))
                lower_lim(i)=floor(min(velP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(velF_kurt(:,i)));
            end

            if(max(velP_kurt(:,i))>max(velF_kurt(:,i)))
                upper_lim(i)=ceil(max(velP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(velF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,velP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,velP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_norm_fig)
        if(min(velP_norm)<min(velF_norm))
            lower_lim(1)=floor(min(velP_norm));
        else
            lower_lim(1)=floor(min(velF_norm));
        end

        if(max(velP_norm)>max(velF_norm))
            upper_lim(1)=ceil(max(velP_norm));
        else
            upper_lim(1)=ceil(max(velF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,velP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,velP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,velF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(velXY_norm_fig)
        if(min(velP_XY_norm)<min(velF_XY_norm))
            lower_lim(1)=floor(min(velP_XY_norm));
        else
            lower_lim(1)=floor(min(velF_XY_norm));
        end

        if(max(velP_XY_norm)>max(velF_XY_norm))
            upper_lim(1)=ceil(max(velP_XY_norm));
        else
            upper_lim(1)=ceil(max(velF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,velP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,velP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,velF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(vel_max_fig)
        for i=1:3
            if(min(velP_max(:,i))<min(velF_max(:,i)))
                lower_lim(i)=floor(min(velP_max(:,i)));
            else
                lower_lim(i)=floor(min(velF_max(:,i)));
            end

            if(max(velP_max(:,i))>max(velF_max(:,i)))
                upper_lim(i)=ceil(max(velP_max(:,i)));
            else
                upper_lim(i)=ceil(max(velF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,velP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,velP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_min_fig)
        for i=1:3
            if(min(velP_min(:,i))<min(velF_min(:,i)))
                lower_lim(i)=floor(min(velP_min(:,i)));
            else
                lower_lim(i)=floor(min(velF_min(:,i)));
            end

            if(max(velP_min(:,i))>max(velF_min(:,i)))
                upper_lim(i)=ceil(max(velP_min(:,i)));
            else
                upper_lim(i)=ceil(max(velF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,velP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,velP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vel_peak_fig)
        for i=1:3
            if(min(velP_peak(:,i))<min(velF_peak(:,i)))
                lower_lim(i)=floor(min(velP_peak(:,i)));
            else
                lower_lim(i)=floor(min(velF_peak(:,i)));
            end

            if(max(velP_peak(:,i))>max(velF_peak(:,i)))
                upper_lim(i)=ceil(max(velP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(velF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,velP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,velP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,velP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,velF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,velP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,velP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,velP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,velF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,velF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,velF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end


% Velocità Angolare
str="Velocità Angolare";
unit="deg/s";

if(vang_fig)
    for i=1:3
        if(min(vangP(:,i))<min(vangF(:,i)))
            lower_lim(i)=floor(min(vangP(:,i)));
        else
            lower_lim(i)=floor(min(vangF(:,i)));
        end

        if(max(vangP(:,i))>max(vangF(:,i)))
            upper_lim(i)=ceil(max(vangP(:,i)));
        else
            upper_lim(i)=ceil(max(vangF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,vangP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,vangP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,vangP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,vangF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,vangF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,vangF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid


    if(vang_mean_fig)
        for i=1:3
            if(min(vangP_mean(:,i))<min(vangF_mean(:,i)))
                lower_lim(i)=floor(min(vangP_mean(:,i)));
            else
                lower_lim(i)=floor(min(vangF_mean(:,i)));
            end

            if(max(vangP_mean(:,i))>max(vangF_mean(:,i)))
                upper_lim(i)=ceil(max(vangP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,vangP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,vangP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_median_fig)
        for i=1:3
            if(min(vangP_median(:,i))<min(vangF_median(:,i)))
                lower_lim(i)=floor(min(vangP_median(:,i)));
            else
                lower_lim(i)=floor(min(vangF_median(:,i)));
            end

            if(max(vangP_median(:,i))>max(vangF_median(:,i)))
                upper_lim(i)=ceil(max(vangP_median(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,vangP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,vangP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(vang_var_fig)
        for i=1:3
            if(min(vangP_var(:,i))<min(vangF_var(:,i)))
                lower_lim(i)=floor(min(vangP_var(:,i)));
            else
                lower_lim(i)=floor(min(vangF_var(:,i)));
            end

            if(max(vangP_var(:,i))>max(vangF_var(:,i)))
                upper_lim(i)=ceil(max(vangP_var(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,vangP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,vangP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_std_fig)
        for i=1:3
            if(min(vangP_std(:,i))<min(vangF_std(:,i)))
                lower_lim(i)=floor(min(vangP_std(:,i)));
            else
                lower_lim(i)=floor(min(vangF_std(:,i)));
            end

            if(max(vangP_std(:,i))>max(vangF_std(:,i)))
                upper_lim(i)=ceil(max(vangP_std(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,vangP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,vangP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_rms_fig)
        for i=1:3
            if(min(vangP_rms(:,i))<min(vangF_rms(:,i)))
                lower_lim(i)=floor(min(vangP_rms(:,i)));
            else
                lower_lim(i)=floor(min(vangF_rms(:,i)));
            end

            if(max(vangP_rms(:,i))>max(vangF_rms(:,i)))
                upper_lim(i)=ceil(max(vangP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,vangP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,vangP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_kurt_fig)
        for i=1:3
            if(min(vangP_kurt(:,i))<min(vangF_kurt(:,i)))
                lower_lim(i)=floor(min(vangP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(vangF_kurt(:,i)));
            end

            if(max(vangP_kurt(:,i))>max(vangF_kurt(:,i)))
                upper_lim(i)=ceil(max(vangP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,vangP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,vangP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_norm_fig)
        if(min(vangP_norm)<min(vangF_norm))
            lower_lim(1)=floor(min(vangP_norm));
        else
            lower_lim(1)=floor(min(vangF_norm));
        end

        if(max(vangP_norm)>max(vangF_norm))
            upper_lim(1)=ceil(max(vangP_norm));
        else
            upper_lim(1)=ceil(max(vangF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,vangP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,vangP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,vangF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(vangXY_norm_fig)
        if(min(vangP_XY_norm)<min(vangF_XY_norm))
            lower_lim(1)=floor(min(vangP_XY_norm));
        else
            lower_lim(1)=floor(min(vangF_XY_norm));
        end

        if(max(vangP_XY_norm)>max(vangF_XY_norm))
            upper_lim(1)=ceil(max(vangP_XY_norm));
        else
            upper_lim(1)=ceil(max(vangF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,vangP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,vangP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,vangF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(vang_max_fig)
        for i=1:3
            if(min(vangP_max(:,i))<min(vangF_max(:,i)))
                lower_lim(i)=floor(min(vangP_max(:,i)));
            else
                lower_lim(i)=floor(min(vangF_max(:,i)));
            end

            if(max(vangP_max(:,i))>max(vangF_max(:,i)))
                upper_lim(i)=ceil(max(vangP_max(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,vangP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,vangP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_min_fig)
        for i=1:3
            if(min(vangP_min(:,i))<min(vangF_min(:,i)))
                lower_lim(i)=floor(min(vangP_min(:,i)));
            else
                lower_lim(i)=floor(min(vangF_min(:,i)));
            end

            if(max(vangP_min(:,i))>max(vangF_min(:,i)))
                upper_lim(i)=ceil(max(vangP_min(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,vangP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,vangP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(vang_peak_fig)
        for i=1:3
            if(min(vangP_peak(:,i))<min(vangF_peak(:,i)))
                lower_lim(i)=floor(min(vangP_peak(:,i)));
            else
                lower_lim(i)=floor(min(vangF_peak(:,i)));
            end

            if(max(vangP_peak(:,i))>max(vangF_peak(:,i)))
                upper_lim(i)=ceil(max(vangP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(vangF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,vangP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,vangP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,vangP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,vangF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,vangP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,vangP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,vangP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,vangF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,vangF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,vangF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel("m/s")
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end

% Angoli
str="Angoli";
unit="rad";

if(ang_fig)
    for i=1:3
        if(min(angP(:,i))<min(angF(:,i)))
            lower_lim(i)=floor(min(angP(:,i)));
        else
            lower_lim(i)=floor(min(angF(:,i)));
        end

        if(max(angP(:,i))>max(angF(:,i)))
            upper_lim(i)=ceil(max(angP(:,i)));
        else
            upper_lim(i)=ceil(max(angF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,angP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,angP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,angP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,angF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,angF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,angF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid


    if(ang_mean_fig)
        for i=1:3
            if(min(angP_mean(:,i))<min(angF_mean(:,i)))
                lower_lim(i)=floor(min(angP_mean(:,i)));
            else
                lower_lim(i)=floor(min(angF_mean(:,i)));
            end

            if(max(angP_mean(:,i))>max(angF_mean(:,i)))
                upper_lim(i)=ceil(max(angP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(angF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,angP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,angP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_median_fig)
        for i=1:3
            if(min(angP_median(:,i))<min(angF_median(:,i)))
                lower_lim(i)=floor(min(angP_median(:,i)));
            else
                lower_lim(i)=floor(min(angF_median(:,i)));
            end

            if(max(angP_median(:,i))>max(angF_median(:,i)))
                upper_lim(i)=ceil(max(angP_median(:,i)));
            else
                upper_lim(i)=ceil(max(angF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,angP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,angP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(ang_var_fig)
        for i=1:3
            if(min(angP_var(:,i))<min(angF_var(:,i)))
                lower_lim(i)=floor(min(angP_var(:,i)));
            else
                lower_lim(i)=floor(min(angF_var(:,i)));
            end

            if(max(angP_var(:,i))>max(angF_var(:,i)))
                upper_lim(i)=ceil(max(angP_var(:,i)));
            else
                upper_lim(i)=ceil(max(angF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,angP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,angP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_std_fig)
        for i=1:3
            if(min(angP_std(:,i))<min(angF_std(:,i)))
                lower_lim(i)=floor(min(angP_std(:,i)));
            else
                lower_lim(i)=floor(min(angF_std(:,i)));
            end

            if(max(angP_std(:,i))>max(angF_std(:,i)))
                upper_lim(i)=ceil(max(angP_std(:,i)));
            else
                upper_lim(i)=ceil(max(angF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,angP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,angP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_rms_fig)
        for i=1:3
            if(min(angP_rms(:,i))<min(angF_rms(:,i)))
                lower_lim(i)=floor(min(angP_rms(:,i)));
            else
                lower_lim(i)=floor(min(angF_rms(:,i)));
            end

            if(max(angP_rms(:,i))>max(angF_rms(:,i)))
                upper_lim(i)=ceil(max(angP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(angF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,angP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,angP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_kurt_fig)
        for i=1:3
            if(min(angP_kurt(:,i))<min(angF_kurt(:,i)))
                lower_lim(i)=floor(min(angP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(angF_kurt(:,i)));
            end

            if(max(angP_kurt(:,i))>max(angF_kurt(:,i)))
                upper_lim(i)=ceil(max(angP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(angF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,angP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,angP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_norm_fig)
        if(min(angP_norm)<min(angF_norm))
            lower_lim(1)=floor(min(angP_norm));
        else
            lower_lim(1)=floor(min(angF_norm));
        end

        if(max(angP_norm)>max(angF_norm))
            upper_lim(1)=ceil(max(angP_norm));
        else
            upper_lim(1)=ceil(max(angF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,angP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,angP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,angF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(angXY_norm_fig)
        if(min(angP_XY_norm)<min(angF_XY_norm))
            lower_lim(1)=floor(min(angP_XY_norm));
        else
            lower_lim(1)=floor(min(angF_XY_norm));
        end

        if(max(angP_XY_norm)>max(angF_XY_norm))
            upper_lim(1)=ceil(max(angP_XY_norm));
        else
            upper_lim(1)=ceil(max(angF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,angP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,angP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,angF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(ang_max_fig)
        for i=1:3
            if(min(angP_max(:,i))<min(angF_max(:,i)))
                lower_lim(i)=floor(min(angP_max(:,i)));
            else
                lower_lim(i)=floor(min(angF_max(:,i)));
            end

            if(max(angP_max(:,i))>max(angF_max(:,i)))
                upper_lim(i)=ceil(max(angP_max(:,i)));
            else
                upper_lim(i)=ceil(max(angF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,angP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,angP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_min_fig)
        for i=1:3
            if(min(angP_min(:,i))<min(angF_min(:,i)))
                lower_lim(i)=floor(min(angP_min(:,i)));
            else
                lower_lim(i)=floor(min(angF_min(:,i)));
            end

            if(max(angP_min(:,i))>max(angF_min(:,i)))
                upper_lim(i)=ceil(max(angP_min(:,i)));
            else
                upper_lim(i)=ceil(max(angF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,angP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,angP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(ang_peak_fig)
        for i=1:3
            if(min(angP_peak(:,i))<min(angF_peak(:,i)))
                lower_lim(i)=floor(min(angP_peak(:,i)));
            else
                lower_lim(i)=floor(min(angF_peak(:,i)));
            end

            if(max(angP_peak(:,i))>max(angF_peak(:,i)))
                upper_lim(i)=ceil(max(angP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(angF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,angP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,angP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,angP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,angF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,angP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,angP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,angP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,angF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,angF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,angF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end


% Campo Magnetico
str="Campo Magnetico";
unit="µT";

if(mag_fig)
    for i=1:3
        if(min(magP(:,i))<min(magF(:,i)))
            lower_lim(i)=floor(min(magP(:,i)));
        else
            lower_lim(i)=floor(min(magF(:,i)));
        end

        if(max(magP(:,i))>max(magF(:,i)))
            upper_lim(i)=ceil(max(magP(:,i)));
        else
            upper_lim(i)=ceil(max(magF(:,i)));
        end
    end

    figure("Name",str)
    subplot(3,2,1)
    plot(tP,magP(:,1),LineWidth=1,Color="r")
    title(str+" Piano")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,3)
    plot(tP,magP(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,5)
    plot(tP,magP(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid

    subplot(3,2,2)
    plot(tF,magF(:,1),LineWidth=1,Color="r")
    title(str+" Forte")
    subtitle("X")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(1),upper_lim(1)])
    grid
    subplot(3,2,4)
    plot(tF,magF(:,2),LineWidth=1,Color="g")
    subtitle("Y")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(2),upper_lim(2)])
    grid
    subplot(3,2,6)
    plot(tF,magF(:,3),LineWidth=1,Color="b")
    subtitle("Z")
    xlabel("t(s)")
    ylabel(unit)
    ylim([lower_lim(3),upper_lim(3)])
    grid


    if(mag_mean_fig)
        for i=1:3
            if(min(magP_mean(:,i))<min(magF_mean(:,i)))
                lower_lim(i)=floor(min(magP_mean(:,i)));
            else
                lower_lim(i)=floor(min(magF_mean(:,i)));
            end

            if(max(magP_mean(:,i))>max(magF_mean(:,i)))
                upper_lim(i)=ceil(max(magP_mean(:,i)));
            else
                upper_lim(i)=ceil(max(magF_mean(:,i)));
            end
        end

        figure("Name",str+" Media")
        subplot(3,1,1)
        plot(tP,magP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_mean(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_mean(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_mean(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_mean(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_mean(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Media")
        subplot(3,2,1)
        plot(tP,magP_mean(:,1),LineWidth=1,Color="r")
        title(str+" Piano Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_mean(:,1),LineWidth=1,Color="r")
        title(str+" Forte Media")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_mean(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_mean(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_median_fig)
        for i=1:3
            if(min(magP_median(:,i))<min(magF_median(:,i)))
                lower_lim(i)=floor(min(magP_median(:,i)));
            else
                lower_lim(i)=floor(min(magF_median(:,i)));
            end

            if(max(magP_median(:,i))>max(magF_median(:,i)))
                upper_lim(i)=ceil(max(magP_median(:,i)));
            else
                upper_lim(i)=ceil(max(magF_median(:,i)));
            end
        end

        figure("Name",str+" Mediana")
        subplot(3,1,1)
        plot(tP,magP_median(:,1),LineWidth=1,Color="r")
        title(str+" Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_median(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_median(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_median(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_median(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_median(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Mediana")
        subplot(3,2,1)
        plot(tP,magP_median(:,1),LineWidth=1,Color="r")
        title(str+" Piano Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_median(:,1),LineWidth=1,Color="r")
        title(str+" Forte Mediana")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_median(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_median(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end


    if(mag_var_fig)
        for i=1:3
            if(min(magP_var(:,i))<min(magF_var(:,i)))
                lower_lim(i)=floor(min(magP_var(:,i)));
            else
                lower_lim(i)=floor(min(magF_var(:,i)));
            end

            if(max(magP_var(:,i))>max(magF_var(:,i)))
                upper_lim(i)=ceil(max(magP_var(:,i)));
            else
                upper_lim(i)=ceil(max(magF_var(:,i)));
            end
        end

        figure("Name",str+" Varianza")
        subplot(3,1,1)
        plot(tP,magP_var(:,1),LineWidth=1,Color="r")
        title(str+" Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_var(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_var(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_var(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_var(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_var(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Varianza")
        subplot(3,2,1)
        plot(tP,magP_var(:,1),LineWidth=1,Color="r")
        title(str+" Piano Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_var(:,1),LineWidth=1,Color="r")
        title(str+" Forte Varianza")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_var(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_var(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_std_fig)
        for i=1:3
            if(min(magP_std(:,i))<min(magF_std(:,i)))
                lower_lim(i)=floor(min(magP_std(:,i)));
            else
                lower_lim(i)=floor(min(magF_std(:,i)));
            end

            if(max(magP_std(:,i))>max(magF_std(:,i)))
                upper_lim(i)=ceil(max(magP_std(:,i)));
            else
                upper_lim(i)=ceil(max(magF_std(:,i)));
            end
        end

        figure("Name",str+" Deviazione Standard")
        subplot(3,1,1)
        plot(tP,magP_std(:,1),LineWidth=1,Color="r")
        title(str+" Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_std(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_std(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_std(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_std(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_std(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Deviazione Standard")
        subplot(3,2,1)
        plot(tP,magP_std(:,1),LineWidth=1,Color="r")
        title(str+" Piano Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_std(:,1),LineWidth=1,Color="r")
        title(str+" Forte Deviazione Standard")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_std(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_std(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_rms_fig)
        for i=1:3
            if(min(magP_rms(:,i))<min(magF_rms(:,i)))
                lower_lim(i)=floor(min(magP_rms(:,i)));
            else
                lower_lim(i)=floor(min(magF_rms(:,i)));
            end

            if(max(magP_rms(:,i))>max(magF_rms(:,i)))
                upper_lim(i)=ceil(max(magP_rms(:,i)));
            else
                upper_lim(i)=ceil(max(magF_rms(:,i)));
            end
        end

        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,1,1)
        plot(tP,magP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_rms(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_rms(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_rms(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_rms(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_rms(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Scarto Quadratico Medio")
        subplot(3,2,1)
        plot(tP,magP_rms(:,1),LineWidth=1,Color="r")
        title(str+" Piano Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_rms(:,1),LineWidth=1,Color="r")
        title(str+" Forte Scarto Quadratico Medio")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_rms(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_rms(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_kurt_fig)
        for i=1:3
            if(min(magP_kurt(:,i))<min(magF_kurt(:,i)))
                lower_lim(i)=floor(min(magP_kurt(:,i)));
            else
                lower_lim(i)=floor(min(magF_kurt(:,i)));
            end

            if(max(magP_kurt(:,i))>max(magF_kurt(:,i)))
                upper_lim(i)=ceil(max(magP_kurt(:,i)));
            else
                upper_lim(i)=ceil(max(magF_kurt(:,i)));
            end
        end

        figure("Name",str+" Kurtosi")
        subplot(3,1,1)
        plot(tP,magP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_kurt(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_kurt(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_kurt(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_kurt(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_kurt(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Kurtosi")
        subplot(3,2,1)
        plot(tP,magP_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Piano Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_kurt(:,1),LineWidth=1,Color="r")
        title(str+" Forte Kurtosi")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_kurt(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_kurt(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_norm_fig)
        if(min(magP_norm)<min(magF_norm))
            lower_lim(1)=floor(min(magP_norm));
        else
            lower_lim(1)=floor(min(magF_norm));
        end

        if(max(magP_norm)>max(magF_norm))
            upper_lim(1)=ceil(max(magP_norm));
        else
            upper_lim(1)=ceil(max(magF_norm));
        end

        figure("Name","Norma "+str)
        plot(tP,magP_norm,LineWidth=1,Color="r")
        title("Norma "+str)
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma "+str)
        subplot(1,2,1)
        plot(tP,magP_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,magF_norm,LineWidth=1,Color="r")
        title("Norma "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(magXY_norm_fig)
        if(min(magP_XY_norm)<min(magF_XY_norm))
            lower_lim(1)=floor(min(magP_XY_norm));
        else
            lower_lim(1)=floor(min(magF_XY_norm));
        end

        if(max(magP_XY_norm)>max(magF_XY_norm))
            upper_lim(1)=ceil(max(magP_XY_norm));
        else
            upper_lim(1)=ceil(max(magF_XY_norm));
        end

        figure("Name","Norma  "+str+" XY")
        plot(tP,magP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_XY_norm,LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")

        figure("Name","Norma  "+str+" XY")
        subplot(1,2,1)
        plot(tP,magP_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

        subplot(1,2,2)
        plot(tF,magF_XY_norm,LineWidth=1,Color="r")
        title("Norma  "+str+" XY Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid

    end

    if(mag_max_fig)
        for i=1:3
            if(min(magP_max(:,i))<min(magF_max(:,i)))
                lower_lim(i)=floor(min(magP_max(:,i)));
            else
                lower_lim(i)=floor(min(magF_max(:,i)));
            end

            if(max(magP_max(:,i))>max(magF_max(:,i)))
                upper_lim(i)=ceil(max(magP_max(:,i)));
            else
                upper_lim(i)=ceil(max(magF_max(:,i)));
            end
        end

        figure("Name",str+" Max")
        subplot(3,1,1)
        plot(tP,magP_max(:,1),LineWidth=1,Color="r")
        title(str+" Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_max(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_max(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_max(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_max(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_max(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Max")
        subplot(3,2,1)
        plot(tP,magP_max(:,1),LineWidth=1,Color="r")
        title(str+" Piano Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_max(:,1),LineWidth=1,Color="r")
        title(str+" Forte Max")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_max(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_max(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_min_fig)
        for i=1:3
            if(min(magP_min(:,i))<min(magF_min(:,i)))
                lower_lim(i)=floor(min(magP_min(:,i)));
            else
                lower_lim(i)=floor(min(magF_min(:,i)));
            end

            if(max(magP_min(:,i))>max(magF_min(:,i)))
                upper_lim(i)=ceil(max(magP_min(:,i)));
            else
                upper_lim(i)=ceil(max(magF_min(:,i)));
            end
        end

        figure("Name",str+" Min")
        subplot(3,1,1)
        plot(tP,magP_min(:,1),LineWidth=1,Color="r")
        title(str+" Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_min(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_min(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_min(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_min(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_min(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name",str+" Min")
        subplot(3,2,1)
        plot(tP,magP_min(:,1),LineWidth=1,Color="r")
        title(str+" Piano Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_min(:,1),LineWidth=1,Color="r")
        title(str+" Forte Min")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_min(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_min(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

    if(mag_peak_fig)
        for i=1:3
            if(min(magP_peak(:,i))<min(magF_peak(:,i)))
                lower_lim(i)=floor(min(magP_peak(:,i)));
            else
                lower_lim(i)=floor(min(magF_peak(:,i)));
            end

            if(max(magP_peak(:,i))>max(magF_peak(:,i)))
                upper_lim(i)=ceil(max(magP_peak(:,i)));
            else
                upper_lim(i)=ceil(max(magF_peak(:,i)));
            end
        end

        figure("Name","Distanza Picco Picco "+str)
        subplot(3,1,1)
        plot(tP,magP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str)
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_peak(:,1),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,2)
        plot(tP,magP_peak(:,2),LineWidth=1,Color="r")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_peak(:,2),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")
        subplot(3,1,3)
        plot(tP,magP_peak(:,3),LineWidth=1,Color="r")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        grid
        hold on
        plot(tF,magF_peak(:,3),LineWidth=1,Color="b")
        legend(str+" Piano", str+" Forte")


        figure("Name","Distanza Picco Picco "+str)
        subplot(3,2,1)
        plot(tP,magP_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Piano")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,3)
        plot(tP,magP_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,5)
        plot(tP,magP_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

        subplot(3,2,2)
        plot(tF,magF_peak(:,1),LineWidth=1,Color="r")
        title("Distanza Picco Picco "+str+" Forte")
        subtitle("X")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(1),upper_lim(1)])
        grid
        subplot(3,2,4)
        plot(tF,magF_peak(:,2),LineWidth=1,Color="g")
        subtitle("Y")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(2),upper_lim(2)])
        grid
        subplot(3,2,6)
        plot(tF,magF_peak(:,3),LineWidth=1,Color="b")
        subtitle("Z")
        xlabel("t(s)")
        ylabel(unit)
        ylim([lower_lim(3),upper_lim(3)])
        grid

    end

end





%% Funzioni

function[sqm] = movrms(f,n)
sqm=zeros(length(f),3);

for i=1:n
    sqm(i,1)=rms([zeros(n-i),f(1:i,1)]);
    sqm(i,2)=rms([zeros(n-i),f(1:i,2)]);
    sqm(i,3)=rms([zeros(n-i),f(1:i,3)]);
end

for i=n+1:length(f)
    sqm(i,1)=rms(f(i-n:i,1));
    sqm(i,2)=rms(f(i-n:i,2));
    sqm(i,3)=rms(f(i-n:i,3));
end

end

function[kurt] = movkurt(f,n)
kurt=zeros(length(f),3);

for i=1:n
    kurt(i,1)=kurtosis([zeros(n-i),f(1:i,1)]);
    kurt(i,2)=kurtosis([zeros(n-i),f(1:i,2)]);
    kurt(i,3)=kurtosis([zeros(n-i),f(1:i,3)]);
end

for i=n+1:length(f)
    kurt(i,1)=kurtosis(f(i-n:i,1));
    kurt(i,2)=kurtosis(f(i-n:i,2));
    kurt(i,3)=kurtosis(f(i-n:i,3));
end

end