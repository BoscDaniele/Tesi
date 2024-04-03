function [t,acc,vang,mag] = AggiustaFrequenza(t_rotta,acc_rotta,vang_rotta,mag_rotta)
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
            disp("lammerda, abbiamo saltato 0.004s "+num2str(i));
            brake;
        end
    end
   
end
end