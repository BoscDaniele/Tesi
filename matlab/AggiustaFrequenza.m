function [t,acc,vang,mag] = AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta)
% Funzione nata dalla necessità di selezionare di quele dei thread che partono nel sensore
% volgiamo tenere i dati. Il sensore infatti faceva saltuariamente partire più thread che 
% che raccoglievano dati ogni 0.04s salvandoli tutti sulla scheda sd. Questo rendeva il tempo di
% acquisizione non più affidabile.

% In questo script vengono prodotti i nuovi vettori t, acc, vang e mag a partire da quelli vecchi
% selezionando i dati che si presentano ogni 0.04s a partire dal primo dato registrato. 

t=[t_rotta(1)];
last=1;

acc=[acc_rotta(1,:)];
vang=[vang_rotta(1,:)];
mag=[mag_rotta(1,:)];

for i=2:length(t_rotta)
    if(t_rotta(i)-t(end)>=0.039 && t_rotta(i)-t(end)<0.041)
        t=[t,t_rotta(i)];
        acc=[acc;acc_rotta(i,:)];
        vang=[vang;vang_rotta(i,:)];
        mag=[mag;mag_rotta(i,:)];
        last=i;

        if(t_rotta(i)-t(end)>0.07)
            disp("abbiamo saltato 0.004s "+num2str(i));
            brake;
        end
    end
   
end
end