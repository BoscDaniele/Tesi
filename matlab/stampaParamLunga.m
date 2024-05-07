clear all
close all
clc

path=".\dati\lunga_forte\";
sr=25; % sample rate

rilievo=2;

[gzRot,gMedio] = GZRot(path);

db=importdata(path + "BlueCoin_Log_N00"+rilievo+".csv").data;

t_rotta=db(:,1)*1e-3;
t_rotta=t_rotta-t_rotta(1);

acc_rotta=db(:,2:4)*gzRot*9.81/-gMedio;
vang_rotta=(db(:,5:7)*1e-3);
mag_rotta=([db(:,8),-db(:,9),db(:,10)]*1e-1);

[t,acc,vang,mag]=AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta);

fun=[{acc(:,1:2)},{vang}];
fun_str=["Acc","V_Ang"];
fun_axes=[{['X','Y']},{['R','P','Y']}];
fun_units=['m/s^2','rad/s'];

for f=1:length(fun)
    funzione=cell2mat(fun(f));

    %% Parametro
    % stampa(t,funzione,fun_str(f),cell2mat(fun_axes(f)),fun_str(f),'t(s)',fun_units(f))

    % %% LowPass
    % fun_low=lowpass(funzione,0.5,25);
    % stampa(t,fun_low,fun_str(f),cell2mat(fun_axes(f)),"LowPass",'t(s)',fun_units(f))


    % %% Media
    % fun_media=movmean(funzione,40);
    % stampa(t,fun_media,fun_str(f),cell2mat(fun_axes(f)),"Media",'t(s)',fun_units(f))


    % %% Valore Medio Rettificato
    % fun_arv=movmean(abs(funzione),40);
    % stampa(t,fun_arv,fun_str(f),cell2mat(fun_axes(f)),"Media Rettificata",'t(s)',fun_units(f))


    % %% Varianza
    % fun_var=movvar(funzione,40);
    % stampa(t,fun_var,fun_str(f),cell2mat(fun_axes(f)),"Varianza",'t(s)',fun_units(f))


    % %% Deviazione Standard
    % fun_std=movstd(funzione,40);
    % stampa(t,fun_std,fun_str(f),cell2mat(fun_axes(f)),"Deviazione Standard",'t(s)',fun_units(f))


    % %% Scarto Quadratico Medio
    % fun_rms=movrms(funzione,40);
    % stampa(t,fun_rms,fun_str(f),cell2mat(fun_axes(f)),"Scarto Quadratico Medio",'t(s)',fun_units(f))


    % %% Kurtosi
    % fun_krt=movkurt(funzione,40);
    % stampa(t,fun_krt,fun_str(f),cell2mat(fun_axes(f)),"Kurtosi",'t(s)',fun_units(f))


    % %% Skewness
    % fun_skw=movskw(funzione,40);
    % stampa(t,fun_skw,fun_str(f),cell2mat(fun_axes(f)),"Skewness",'t(s)',fun_units(f))


    % %% Shape Factor
    % fun_rms=movrms(funzione,40);
    % fun_arv=movmean(abs(funzione),40);
    %
    % fun_shf=fun_rms./fun_arv;
    % stampa(t,fun_shf,fun_str(f),cell2mat(fun_axes(f)),"Shape Factor",'t(s)',fun_units(f))


    % %% Crest Factor
    % fun_max=movmax(funzione,40);
    % fun_rms=movrms(funzione,40);
    %
    % fun_crf=fun_max./fun_rms;
    % stampa(t,fun_crf,fun_str(f),cell2mat(fun_axes(f)),"Crest Factor",'t(s)',fun_units(f))


    % %% Impulse Factor
    % fun_max=movmax(funzione,40);
    % fun_arv=movmean(abs(funzione),40);
    %
    % fun_impf=fun_max./fun_arv;
    % stampa(t,fun_impf,fun_str(f),cell2mat(fun_axes(f)),"Impulse Factor",'t(s)',fun_units(f))


    % %% Margin Factor
    % fun_max=movmax(funzione,40);
    % fun_arvq=movmean(sqrt(abs(funzione)),40).^2;
    %
    % fun_mrgf=fun_max./fun_arvq;
    % stampa(t,fun_mrgf,fun_str(f),cell2mat(fun_axes(f)),"Margin Factor",'t(s)',fun_units(f))


    % %% Max
    % n_max=10;
    % fun_max=movmax(funzione,n_max);
    % stampa(t,fun_max,fun_str(f),cell2mat(fun_axes(f)),"Max",'t(s)',fun_units(f))


    % %% Min
    % n_min=10;
    % fun_min=movmin(funzione,n_min);
    %
    % n=length(cell2mat(fun_axes(f)));
    % f_min=floor(min(funzione));
    % f_max=ceil(max(funzione));
    % limY=[f_min',f_max'];
    %
    % textLimY=ones(n,5);
    %
    % for i=1:n
    %     textLimY(i,:)=(limY(i,2))*ones(5,1);
    % end
    %
    % limY(:,2)=(limY(:,2)+1).*1.25;
    %
    % stampa_gen(t,fun_min,fun_str(f),cell2mat(fun_axes(f)),"Min",'t(s)',fun_units(f),limY,textLimY)


    % %% Peak
    % n_Peak=10;
    % fun_max=movmax(funzione,n_Peak);
    % fun_min=movmin(funzione,n_Peak);
    %
    % fun_peak=fun_max-fun_min;
    % stampa(t,fun_peak,fun_str(f),cell2mat(fun_axes(f)),"Peak",'t(s)',fun_units(f))


    %% Trasformata
    sezione=ones(5,2);
    sezione(1,:)=[1 15];
    sezione(2,:)=[15 19.4];
    sezione(3,:)=[19.5 23.76];
    sezione(4,:)=[23.76 25.4];
    sezione(5,:)=[25.4 29.76];

    str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"];

    L=length(funzione);
    frequenza=sr/L*(0:(L/2));

    Y_noMedia=fft(funzione-movmean(funzione,40));
    P2_noMedia=abs(Y_noMedia/L);
    trasformata=P2_noMedia(1:(L/2+1),:);
    trasformata(2:end-1,:)=2*trasformata(2:end-1,:);

    % stampa_freq(frequenza,trasformata,cell2mat(fun_axes(f)),"Trasformata")


    xdftP=Y_noMedia(1:L/2+1,:);
    spettro=(1/(sr*L))*abs(xdftP).^2;
    spettro(2:end-1,:)=2*spettro(2:end-1,:);

    stampa_freq(frequenza,spettro,cell2mat(fun_axes(f)),"Spettro")


    trasform_mean=zeros(length(sezione),length(cell2mat(fun_axes(f))));
    trasform_cent=zeros(length(sezione),length(cell2mat(fun_axes(f))));
    trasform_var=zeros(length(sezione),length(cell2mat(fun_axes(f))));
    entropia=zeros(length(sezione),length(cell2mat(fun_axes(f))));

    for i=1:length(sezione)
        L=(sezione(i,2)-sezione(i,1))*25;
        frequenza=sr/L*(0:(L/2));

        Y_noMedia=fft(funzione(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(funzione(sr*sezione(i,1):sr*sezione(i,2),:),40));
        P2_noMedia=abs(Y_noMedia/L);
        trasformata=P2_noMedia(1:(L/2+1),:);
        trasformata(2:end-1,:)=2*trasformata(2:end-1,:);

        trasform_mean(i,:)=mean(trasformata);
        trasform_cent(i,:)=mean(trasformata.*frequenza');

        for j=1:length(cell2mat(fun_axes(f)))
            trasform_var(i,j)=sum((frequenza'-trasform_cent(i,j)).*trasformata(:,j))/sum(trasformata(:,j));
        end

        % stampa_freq(frequenza,trasformata,cell2mat(fun_axes(f)),"Trasformata "+str(i))


        xdftP=Y_noMedia(1:L/2+1,:);
        spettro=(1/(sr*L))*abs(xdftP).^2;
        spettro(2:end-1,:)=2*spettro(2:end-1,:);

        % stampa_freq(frequenza,spettro,cell2mat(fun_axes(f)),"Spettro "+str(i))

        for j=1:length(cell2mat(fun_axes(f)))
            p=spettro(:,j)/sum(spettro(:,j));
            entropia(i,j)=-sum(p.*log2(p));
        end

    end

    % %% Ampiezza Media
    % stampa_freqAmp(trasform_mean,cell2mat(fun_axes(f)),"Ampiezza Media",str)


    % %% Frequency Centroid
    % stampa_freqParam(trasform_cent,cell2mat(fun_axes(f)),"Frequency Centroid",str)


    % %% Frequency Variance
    % stampa_freqParam(trasform_var,cell2mat(fun_axes(f)),"Frequency Variance",str)


    % %% Spectral Entropy
    % stampa_freqParam(entropia,cell2mat(fun_axes(f)),"Spectral Entropy",str)

end




%% Funzioni
function[sqm] = movrms(f,n)
w=width(f);
newf=[zeros(n/2,w);f;zeros(n/2,w)];

sqm=zeros(length(f),w);

for i=1:length(f)
    sqm(i,:)=rms(newf(i:i+n,:));
end

end

function[kurt] = movkurt(f,n)
w=width(f);
newf=[zeros(n/2,w);f;zeros(n/2,w)];

kurt=zeros(length(f),w);

for i=1:length(f)
    kurt(i,:)=kurtosis(newf(i:i+n,:));
end

end

function[skew] = movskw(f,n)
w=width(f);
newf=[zeros(n/2,w);f;zeros(n/2,w)];

skew=zeros(length(f),w);

for i=1:length(f)
    skew(i,:)=skewness(newf(i:i+n,:));
end

end


%% Grafici
function stampa(t,fun,fun_str,fun_axes,tit,xlbl,ylbl)
n=length(fun_axes);

fun_min=floor(min(fun));
fun_max=ceil(max(fun));
limY=[fun_min',fun_max'];

textLimY=ones(n,5);

for i=1:n
    textLimY(i,:)=(limY(i,2)+1)*ones(5,1);
end

limY(:,2)=(limY(:,2)+1).*1.25;

stampa_gen(t,fun,fun_str,fun_axes,tit,xlbl,ylbl,limY,textLimY)
end

function stampa_gen(t,fun,fun_str,fun_axes,tit,xlbl,ylbl,limY,textLimY)
n=length(fun_axes);
font="Times New Roman";

linee=ones(5,2);
linee(1,:)=[1 1];
linee(2,:)=[15 15];
linee(3,:)=[19.4 19.4];
linee(4,:)=[23.76 23.76];
linee(5,:)=[25.4 25.4];
% linee(6,:)=[29.76 29.76];

str={'Accelerate','Idle','Accelerate','Idle','Brake + Tourning'};

for i=1:n
    if (fun_axes(i)=='X')
        c(i)='r';
    elseif (fun_axes(i)=='Y')
        c(i)='g';
    else
        c(i)='b';
    end
end

f=figure;
for i=1:n
    subplot(n,1,i)
    plot(t,fun(:,i),LineWidth=1,color=c(i))
    if(i==1)
        title(tit,FontName=font)
    end
    subtitle(fun_axes(i),FontName=font)
    xlabel(xlbl,FontName=font)
    ylabel(ylbl,FontName=font)
    ylim(limY(i,:))
    grid
    hold on
    for j=1:length(linee)
        line(linee(j,:),limY(i,:), 'Color','black')
    end
    text(linee(:,1)+0.5*[1,1,1,1/2,1]',textLimY(i,:),str,FontSize=6.5, FontWeight="bold",FontName=font)

end

% exportgraphics(f,"D:\Users\Daniele\Desktop\"+fun_str+"\"+tit+".png")

end

function stampa_freq(freq,trasform,fun_axes,tit)
font="Times New Roman";

f=figure;
for i=1:length(fun_axes)
    subplot(length(fun_axes),1,i)
    plot(freq,trasform(:,i),LineWidth=1,Color="r")
    if(i == 1)
        title(tit,FontName=font)
    end
    subtitle(fun_axes(i),FontName=font)
    xlabel('Hz',FontName=font)
    ylabel(fun_axes(i)+"(Hz)",FontName=font)
    grid
end

% exportgraphics(f,"D:\Users\Daniele\Desktop\Trasformata\"+tit+".png")

end

function stampa_freqAmp(param,fun_axes,tit,sezioni)
n=length(sezioni);
font="Times New Roman";

for j=1:length(fun_axes)
    limY=[floor(min(param(:,j))),ceil(max(param(:,j)))];

    f=figure;
    for i=1:n
        plot(1:100,param(i,j)*ones(100,1),LineWidth=1,DisplayName=sezioni(i))
        ylabel(fun_axes(j)+"(Hz)",FontName=font)
        ylim(limY)
        grid
        hold on
    end
    title(tit+" "+fun_axes(j),FontName=font)
    legend
    
    % exportgraphics(f,"D:\Users\Daniele\Desktop\Trasformata\"+tit+" "+fun_axes(j)+".png")
end

end

function stampa_freqParam(param,fun_axes,tit,sezioni)
n=length(sezioni);
font="Times New Roman";

for j=1:length(fun_axes)
    f=figure;
    for i=1:n
        plot(param(i,j)*ones(100,1),0:99,LineWidth=1,DisplayName=sezioni(i))
        ylabel(fun_axes(j)+"(Hz)",FontName=font)
        ylim([0,100])
        grid
        hold on
    end
    title(tit+" "+fun_axes(j),FontName=font)
    legend
    
    % exportgraphics(f,"D:\Users\Daniele\Desktop\Trasformata\"+tit+" "+fun_axes(j)+".png")
end

end




function stampa_SpEnt(fun)
font="Times New Roman";
sr = 25; %sample rate

sezione=ones(5,2);
sezione(1,:)=[1 15];
sezione(2,:)=[15 19.4];
sezione(3,:)=[19.5 23.76];
sezione(4,:)=[23.76 25.4];
sezione(5,:)=[25.4 29.76];

str=["Accelerazione 1","Idle 1","Accelerazione 2","Idle 2","Brake"]';

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    xdftP=Y_noMedia(1:L/2+1,:);
    spettro=(1/(sr*L))*abs(xdftP).^2;
    spettro(2:end-1,:)=2*spettro(2:end-1,:);

    p=zeros(length(spettro),1);
    spettro_tot=sum(spettro(:,1));

    for j=1:length(spettro)
        p(j)=spettro(j,1)/spettro_tot;
    end

    entropia=-sum(p.*log2(p));

    plot(freq,entropia*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("X''(Hz)",FontName=font)
    grid
    hold on
end
title("Spectral Entropy X",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\SpEntX"+str(i)+".png")

f=figure;
for i=1:5
    L=(sezione(i,2)-sezione(i,1))*25;
    freq=sr/L*(0:(L/2));

    Y_noMedia=fft(fun(sr*sezione(i,1):sr*sezione(i,2),:)-movmean(fun(sr*sezione(i,1):sr*sezione(i,2),:),40));
    xdftP=Y_noMedia(1:L/2+1,:);
    spettro=(1/(sr*L))*abs(xdftP).^2;
    spettro(2:end-1,:)=2*spettro(2:end-1,:);

    p=zeros(length(spettro),1);
    spettro_tot=sum(spettro(:,2));

    for j=1:length(spettro)
        p(j)=spettro(j,2)/spettro_tot;
    end

    entropia=-sum(p.*log2(p));

    plot(freq,entropia*ones(length(freq),1),LineWidth=1,DisplayName=str(i))
    xlabel('Hz',FontName=font)
    ylabel("Y''(Hz)",FontName=font)
    grid
    hold on
end
title("Spectral Entropy Y",FontName=font)
legend

% exportgraphics(f,"D:\Users\Daniele\Desktop\SpEntY"+str(i)+".png")
end



