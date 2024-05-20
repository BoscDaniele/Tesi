clear all
close all
clc

pathCurvaU=".\dati\curvaU_forte\";
pathLunga=".\dati\lunga_forte\";

sr=25; % sample rate

rilievoCurvaU=4;
rilievoLunga=2;

%% CurvaU
[gzRotCurvaU,gMedioCurvaU] = GZRot(pathCurvaU);

dbCurvaU=importdata(pathCurvaU + "BlueCoin_Log_N00"+rilievoCurvaU+".csv").data;

tCurvaU_rotta=dbCurvaU(:,1)*1e-3;
tCurvaU_rotta=tCurvaU_rotta-tCurvaU_rotta(1);

accCurvaU_rotta=dbCurvaU(:,2:4)*gzRotCurvaU*9.81/-gMedioCurvaU;
vangCurvaU_rotta=(dbCurvaU(:,5:7)*1e-3);
magCurvaU_rotta=([dbCurvaU(:,8),-dbCurvaU(:,9),dbCurvaU(:,10)]*1e-1)*gzRotCurvaU;

[tCurvaU,accCurvaU,vangCurvaU,magCurvaU]=AggiustaFrequenza(tCurvaU_rotta,accCurvaU_rotta,vangCurvaU_rotta,magCurvaU_rotta);
velCurvaU=cumsum(accCurvaU)*0.04;


%% Lunga
[gzRotLunga,gMedioLunga] = GZRot(pathLunga);

dbLunga=importdata(pathLunga + "BlueCoin_Log_N00"+rilievoLunga+".csv").data;

tLunga_rotta=dbLunga(:,1)*1e-3;
tLunga_rotta=tLunga_rotta-tLunga_rotta(1);

accLunga_rotta=dbLunga(:,2:4)*gzRotLunga*9.81/-gMedioLunga;
vangLunga_rotta=(dbLunga(:,5:7)*1e-3);
magLunga_rotta=([dbLunga(:,8),-dbLunga(:,9),dbLunga(:,10)]*1e-1)*gzRotLunga;

[tLunga,accLunga,vangLunga,magLunga]=AggiustaFrequenza(tLunga_rotta,accLunga_rotta,vangLunga_rotta,magLunga_rotta);
velLunga=cumsum(accLunga)*0.04;

%% Fun
acc={accLunga,accCurvaU};
vang={vangLunga,vangCurvaU};
mag={magLunga,magCurvaU};
vel={velLunga,velCurvaU};

fun={acc,vang,mag,vel};
% fun={acc,vang,mag};
fun_str=["Acc","VAng","Mag","Vel"];
fun_axes={["X","Y","Z"],["Roll","Pitch","Yaw"],["X","Y","Z"],["X","Y","Z"]};
fun_units=["m/s^2","deg/s","µT","m/s"];

for f=1:length(fun)
    funzioneLunga=fun{f}{1};
    funzioneCurvaU=fun{f}{2};
    axes=fun_axes{f};

    %% Parametro
    stampa(tLunga,tCurvaU,funzioneLunga,funzioneCurvaU,fun_str(f),axes,fun_str(f),'t(s)',fun_units(f))


    %% LowPass
    funP_low=lowpass(funzioneLunga,0.5,25);
    funF_low=lowpass(funzioneCurvaU,0.5,25);
    stampa(tLunga,tCurvaU,funP_low,funF_low,fun_str(f),axes,"LowPass",'t(s)',fun_units(f))


    %% Media
    funP_media=movmean(funzioneLunga,40);
    funF_media=movmean(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_media,funF_media,fun_str(f),axes,"Media",'t(s)',fun_units(f))


    %% Valore Medio Rettificato
    funP_arv=movmean(abs(funzioneLunga),40);
    funF_arv=movmean(abs(funzioneCurvaU),40);
    stampa(tLunga,tCurvaU,funP_arv,funF_arv,fun_str(f),axes,"Media Rettificata",'t(s)',fun_units(f))


    %% Varianza
    funP_var=movvar(funzioneLunga,40);
    funF_var=movvar(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_var,funF_var,fun_str(f),axes,"Varianza",'t(s)',fun_units(f))


    %% Deviazione Standard
    funP_std=movstd(funzioneLunga,40);
    funF_std=movstd(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_std,funF_std,fun_str(f),axes,"Deviazione Standard",'t(s)',fun_units(f))


    %% Scarto Quadratico Medio
    funP_rms=movrms(funzioneLunga,40);
    funF_rms=movrms(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_rms,funF_rms,fun_str(f),axes,"Scarto Quadratico Medio",'t(s)',fun_units(f))


    %% Kurtosi
    funP_krt=movkurt(funzioneLunga,40);
    funF_krt=movkurt(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_krt,funF_krt,fun_str(f),axes,"Kurtosi",'t(s)',fun_units(f))


    %% Skewness
    funP_skw=movskw(funzioneLunga,40);
    funF_skw=movskw(funzioneCurvaU,40);
    stampa(tLunga,tCurvaU,funP_skw,funF_skw,fun_str(f),axes,"Skewness",'t(s)',fun_units(f))


    %% Max
    n_max=10;
    funP_max=movmax(funzioneLunga,n_max);
    funF_max=movmax(funzioneCurvaU,n_max);
    stampa(tLunga,tCurvaU,funP_max,funF_max,fun_str(f),axes,"Max",'t(s)',fun_units(f))


    %% Min
    n_min=10;
    funP_min=movmin(funzioneLunga,n_min);
    funF_min=movmin(funzioneCurvaU,n_min);

    stampa(tLunga,tCurvaU,funP_min,funF_min,fun_str(f),axes,"Min",'t(s)',fun_units(f))


    %% Peak
    funP_peak=funP_max-funP_min;
    funF_peak=funF_max-funF_min;
    stampa(tLunga,tCurvaU,funP_peak,funF_peak,fun_str(f),axes,"Peak",'t(s)',fun_units(f))


    %% Shape Factor
    funP_shf=funP_rms./funP_arv;
    funF_shf=funF_rms./funF_arv;
    stampa(tLunga,tCurvaU,funP_shf,funF_shf,fun_str(f),axes,"Shape Factor",'t(s)',fun_units(f))


    %% Crest Factor
    funP_crf=funP_max./funP_rms;
    funF_crf=funF_max./funF_rms;
    stampa(tLunga,tCurvaU,funP_crf,funF_crf,fun_str(f),axes,"Crest Factor",'t(s)',fun_units(f))


    %% Impulse Factor
    funP_impf=funP_max./funP_arv;
    funF_impf=funF_max./funF_arv;
    stampa(tLunga,tCurvaU,funP_impf,funF_impf,fun_str(f),axes,"Impulse Factor",'t(s)',fun_units(f))


    %% Margin Factor
    funP_arvq=movmean(sqrt(abs(funzioneLunga)),40).^2;
    funF_arvq=movmean(sqrt(abs(funzioneCurvaU)),40).^2;

    funP_mrgf=funP_max./funP_arvq;
    funF_mrgf=funF_max./funF_arvq;
    stampa(tLunga,tCurvaU,funP_mrgf,funF_mrgf,fun_str(f),axes,"Margin Factor",'t(s)',fun_units(f))


    %% Trasformata
    LP=length(funzioneLunga);
    frequenzaP=sr/LP*(0:(LP/2));

    LF=length(funzioneCurvaU);
    frequenzaF=sr/LF*(0:(LF/2));

    YP_noMedia=fft(funzioneLunga-movmean(funzioneLunga,40));
    P2P_noMedia=abs(YP_noMedia/LP);
    trasformataP=P2P_noMedia(1:(LP/2+1),:);
    trasformataP(2:end-1,:)=2*trasformataP(2:end-1,:);

    YF_noMedia=fft(funzioneCurvaU-movmean(funzioneCurvaU,40));
    P2F_noMedia=abs(YF_noMedia/LF);
    trasformataF=P2F_noMedia(1:(LF/2+1),:);
    trasformataF(2:end-1,:)=2*trasformataF(2:end-1,:);

    stampa_freq(frequenzaP,frequenzaF,trasformataP,trasformataF,axes,fun_str(f),"Trasformata")


    %% Spettro
    xdftP=YP_noMedia(1:LP/2+1,:);
    spettroP=(1/(sr*LP))*abs(xdftP).^2;
    spettroP(2:end-1,:)=2*spettroP(2:end-1,:);

    xdftF=YF_noMedia(1:LF/2+1,:);
    spettroF=(1/(sr*LF))*abs(xdftF).^2;
    spettroF(2:end-1,:)=2*spettroF(2:end-1,:);

    stampa_freq(frequenzaP,frequenzaF,spettroP,spettroF,axes,fun_str(f),"Spettro")


    %% Ampiezza Media
    trasformP_mean=mean(trasformataP);
    trasformF_mean=mean(trasformataF);

    stampa_freqAmp(trasformP_mean,trasformF_mean,axes,fun_str(f),"Ampiezza Media")


    %% Frequency Centroid
    trasformP_cent=mean(trasformataP.*frequenzaP');
    trasformF_cent=mean(trasformataF.*frequenzaF');

    stampa_freqParam(trasformP_cent,trasformF_cent,axes,fun_str(f),"Frequency Centroid")


    %% Frequency Variance
    trasformP_var=zeros(length(axes));
    trasformF_var=zeros(length(axes));

    for i=1:length(axes)
        trasformP_var(i)=sum((frequenzaP'-trasformP_cent(i)).*trasformataP(:,i))/sum(trasformataP(:,i));
        trasformF_var(i)=sum((frequenzaF'-trasformF_cent(i)).*trasformataF(:,i))/sum(trasformataF(:,i));
    end

    stampa_freqParam(trasformP_var,trasformF_var,axes,fun_str(f),"Frequency Variance")


    %% Spectral Entropy
    entropiaP=zeros(length(axes));
    entropiaF=zeros(length(axes));

    for i=1:length(axes)
        pP=spettroP(:,i)/sum(spettroP(:,i));
        pF=spettroF(:,i)/sum(spettroF(:,i));

        entropiaP(i)=-sum(pP.*log2(pP));
        entropiaF(i)=-sum(pF.*log2(pF));
    end

    stampa_freqParam(entropiaP,entropiaF,axes,fun_str(f),"Spectral Entropy")

    close all
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
function stampa(tP,tF,funP,funF,fun_str,fun_axes,tit,xlbl,ylbl)
n=length(fun_axes);

% fun_min=floor(min([funP;funF]));
% fun_max=ceil(max([funP;funF]));

fun_min=(min([funP;funF]));
fun_max=(max([funP;funF]));

limY=[fun_min',fun_max'];

stampa_gen(tP,tF,funP,funF,fun_str,fun_axes,tit,xlbl,ylbl,limY)
end

function stampa_gen(tP,tF,funP,funF,fun_str,fun_axes,tit,xlbl,ylbl,limY)
n=length(fun_axes);
font="Times New Roman";

for i=1:n
    if (fun_axes(i)=="X"||fun_axes(i)=="Roll")
        c(i)='r';
    elseif (fun_axes(i)=="Y"||fun_axes(i)=="Pitch")
        c(i)='g';
    else
        c(i)='b';
    end
end

f=figure;
for i=3:2:2*n+1
    j=floor(i/2);
    subplot(n,2,i-2)
    plot(tP,funP(:,j),LineWidth=1,color=c(j))
    if(i==3)
        title(tit+" Rettilineo",FontName=font)
    end
    subtitle(fun_axes(j),FontName=font)
    xlabel(xlbl,FontName=font)
    ylabel(ylbl,FontName=font)
    ylim(limY(j,:))
    grid
end

for i=2:2:2*n
    subplot(n,2,i)
    plot(tF,funF(:,i/2),LineWidth=1,color=c(i/2))
    if(i==2)
        title(tit+" Curva U",FontName=font)
    end
    subtitle(fun_axes(i/2),FontName=font)
    xlabel(xlbl,FontName=font)
    ylabel(ylbl,FontName=font)
    ylim(limY(i/2,:))
    grid
end

exportgraphics(f,"..\slide\lunga_curvaU\figure\"+fun_str+"\"+tit+".png")

end

function stampa_freq(freqP,freqF,trasformP,trasformF,fun_axes,fun_str,tit)
font="Times New Roman";
n=length(fun_axes);

% limY_min=floor(min([trasformP;trasformF]));
% limY_max=ceil(max([trasformP;trasformF]));

limY_min=(min([trasformP;trasformF]));
limY_max=(max([trasformP;trasformF]));

limY=[limY_min',limY_max'];

f=figure;
for i=3:2:2*n+1
    j=floor(i/2);
    subplot(n,2,i-2)
    p=plot(freqP,trasformP(:,j),LineWidth=1);
    if(j == 1)
        title(tit+" Rettilineo",FontName=font)
    end
    subtitle(fun_axes(j),FontName=font)
    xlabel('Hz',FontName=font)
    ylabel(fun_axes(j)+"(Hz)",FontName=font)
    ylim(limY(j,:))
    grid

    if (fun_axes(j)=="X"||fun_axes(j)=="Roll")
        p.Color='r';
    elseif (fun_axes(j)=="Y"||fun_axes(j)=="Pitch")
        p.Color='g';
    else
        p.Color='b';
    end
end

for i=2:2:2*n
    subplot(n,2,i)
    p=plot(freqF,trasformF(:,i/2),LineWidth=1);
    if(i == 2)
        title(tit+"Curva U",FontName=font)
    end
    subtitle(fun_axes(i/2),FontName=font)
    xlabel('Hz',FontName=font)
    ylabel(fun_axes(i/2)+"(Hz)",FontName=font)
    ylim(limY(i/2,:))
    grid

    if (fun_axes(i/2)=="X"||fun_axes(i/2)=="Roll")
        p.Color='r';
    elseif (fun_axes(i/2)=="Y"||fun_axes(i/2)=="Pitch")
        p.Color='g';
    else
        p.Color='b';
    end
end

exportgraphics(f,"..\slide\lunga_curvaU\figure\"+fun_str+"\Trasformata\"+tit+".png")

end

function stampa_freqAmp(paramP,paramF,fun_axes,fun_str,tit)
font="Times New Roman";

for i=1:length(fun_axes)
    f=figure;
    plot(1:100,paramP(i)*ones(100,1),LineWidth=1,DisplayName=tit+" Rettilineo");
    title(tit+" "+fun_axes(i),FontName=font)
    ylabel(fun_axes(i)+"(Hz)",FontName=font)
    grid
    hold on
    plot(1:100,paramF(i)*ones(100,1),LineWidth=1,DisplayName=tit+" Curva U");

    legend
    exportgraphics(f,"..\slide\lunga_curvaU\figure\"+fun_str+"\Trasformata\"+tit+fun_axes(i)+".png")
end

end

function stampa_freqParam(paramP,paramF,fun_axes,fun_str,tit)
font="Times New Roman";

for i=1:length(fun_axes)
    f=figure;
    plot(paramP(i)*ones(100,1),0:99,LineWidth=1,DisplayName=tit+" Rettilineo")
    title(tit+" "+fun_axes(i),FontName=font)
    ylabel(fun_axes(i)+"(Hz)",FontName=font)
    ylim([0,100])
    grid
    hold on
    plot(paramF(i)*ones(100,1),0:99,LineWidth=1,DisplayName=tit+" Curva U")

    legend
    exportgraphics(f,"..\slide\lunga_curvaU\figure\"+fun_str+"\Trasformata\"+tit+fun_axes(i)+".png")

end

end


