clear
close all
clc

%Questo script serve per controllare che il tempo di campionamento dei
%rilevamenti fatti durante un'uscita sia effettivamente di 0.04s, cosa, a
%quanto pare, non scontata

path="..\db\secchia\";
nRilievi=8;

for i=0:nRilievi-1
    db=importdata(path + "BlueCoin_Log_N00"+i+".csv").data;

    t=db(1:length(db),1)*1e-3;
    t=t-t(1);

    for j=2:length(t)
    intervalloT(j)=t(j)-t(j-1);
    end
    disp("");
    disp("---------------------------------------------");
    disp("tempo di campionamento minimo: "+num2str(min(intervalloT(2:end))));
    disp("tempo di campionamento massimo: "+num2str(max(intervalloT(2:end))));
end

