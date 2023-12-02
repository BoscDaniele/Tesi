close all
clear
clc

%import dati accelerazioni
dbA=importdata("dbdm\BlueCoin_Log_N000.csv");
acc=dbA.data(2:end,2:4);

%estrazione tempi
t=dbA.data(2:end,1)*1e-3;
t=t-t(1);

% %plot
% figure(Name="Accelereazione")
% subplot(3,1,1);
% plot(t,acc(:,1),LineWidth=1,Color="b");
% subtitle("X", Color="b");
% subplot(3,1,2);
% plot(t,acc(:,2),LineWidth=1,Color="r");
% subtitle("Y", Color="b");
% subplot(3,1,3);
% plot(t,acc(:,3),LineWidth=1,Color="g");
% subtitle("Z", Color="b");

% %grafico 3D vettore gravità
% figure
% plot3([0,acc(1,1)],[0,acc(1,2)],[0,acc(1,3)],LineWidth=1,Color="black")
% grid
% hold on
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="b")
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="r")
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="g")
% xlabel("X", Color="b")
% ylabel("Y",Color="r")
% zlabel("Z", Color="g")

% figure(Name="Acc XZ Cavalletto")
% plot([0,acc(1,1)],[0,acc(1,3)]);
% grid;
% xlabel("X", Color="b")
% ylabel("Z",Color="g")

modXZ=sqrt(acc(1,1)^2+acc(1,3)^2);
thetaXZ=-acos(-acc(1,3)/modXZ);
mRotCav=[cos(thetaXZ), 0, sin(thetaXZ); 0, 1, 0; -sin(thetaXZ), 0, cos(thetaXZ)];
accNoCav=acc(1,:)*mRotCav;

% figure(Name="Acc XZ No Cavalletto")
% plot([0,accNoCav(1,1)],[0,accNoCav(1,3)]);
% grid
% xlabel("X", Color="b")
% ylabel("Z",Color="g")

% figure(Name="Acc No Cavalletto")
% subplot(1,3,1);
% grid;
% plot([0,accNoCav(1,1)],[0,accNoCav(1,2)],LineWidth=1);
% subtitle("Acc XY")
% subplot(1,3,2);
% grid;
% plot([0,accNoCav(1,2)],[0,accNoCav(1,3)],LineWidth=1);
% subtitle("Acc YZ")
% subplot(1,3,3);
% grid;
% plot([0,accNoCav(1,1)],[0,accNoCav(1,3)],LineWidth=1);
% subtitle("Acc XZ")

newModYZ=sqrt(accNoCav(1,2)^2+accNoCav(1,3)^2);
newThetaYZ=-acos(-accNoCav(1,3)/newModYZ);
mRot=[ 1, 0, 0; 0, cos(newThetaYZ), -sin(newThetaYZ); 0, sin(newThetaYZ), cos(newThetaYZ)];
mRot=mRotCav*mRot*[cos(pi), -sin(pi), 0; sin(pi), cos(pi), 0; 0, 0, 1];
% mRot=mRotCav;
newAccNoCav=accNoCav(1,:)*mRot;

% figure(Name="New Acc No Cav")
% subplot(1,3,1);
% grid;
% plot([0,newAccNoCav(1,1)],[0,newAccNoCav(1,2)],LineWidth=1);
% subtitle("Acc XY")
% subplot(1,3,2);
% grid;
% plot([0,newAccNoCav(1,2)],[0,newAccNoCav(1,3)],LineWidth=1);
% subtitle("Acc YZ")
% subplot(1,3,3);
% grid;
% plot([0,newAccNoCav(1,1)],[0,newAccNoCav(1,3)],LineWidth=1);
% subtitle("Acc XZ")

%grafico 3D vettore gravità
% figure
% plot3([0,newAccNoCav(1,1)],[0,newAccNoCav(1,2)],[0,newAccNoCav(1,3)],LineWidth=1,Color="black")
% grid
% hold on
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="b")
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="r")
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="g")
% xlabel("X", Color="b")
% ylabel("Y",Color="r")
% zlabel("Z", Color="g")

%import dati accelerazioni
dbA=importdata("dbdm\BlueCoin_Log_N004.csv");
acc=dbA.data(2:end,2:4);

%estrazione tempi
t=dbA.data(2:end,1)*1e-3;
t=t-t(1);

%calcolo nuova matrice delle accelerazioni
newAcc = zeros(length(acc),3);
for i=1:length(acc)
    newAcc(i,:)=acc(i,:)*mRot;
end


% angY=-real(acos(acc(:,3)/-1000));

% for i=(length(angY)+1)*2/3:length(angY)
%     angY(i)=-angY(i);
% end

% disp("angolo rotazione medio: " + num2str(mean(angY)*360/(2*pi)));
% 
% for i=1:length(acc)
%     acc(i,:)=acc(i,:)*[cos(angY(i)), 0, sin(angY(i)); 0, 1, 0; -sin(angY(i)), 0, cos(angY(i))];
% end
% newAcc(:,1)=newAcc(:,1)+xcomponent;

acc=acc.*(9.81/1000);

newAcc=newAcc.*(9.81/1000);


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

%plot
figure(Name="Accelereazione")
subplot(3,1,1);
plot(t,newAcc(:,1),LineWidth=1,Color="b");
subtitle("X", Color="b");
subplot(3,1,2);
plot(t,newAcc(:,2),LineWidth=1,Color="r");
subtitle("Y", Color="b");
subplot(3,1,3);
plot(t,newAcc(:,3),LineWidth=1,Color="g");
subtitle("Z", Color="b");

% %grafico accellerazione in x
% figure
% plot(t,acc(:,1))
% grid on
% title('Accelerazione in X')

%% Cacolo velocità

%velocità in x
vx=cumtrapz(t,acc(:,1));
newVx=cumtrapz(t,newAcc(:,1));

%grafico velocità in x
figure(Name="Asse X")
subplot(1,2,1)
plot(t,acc(:,1))
grid on
subtitle('Accelerazione')
subplot(1,2,2)
plot(t,vx(:))
grid on
subtitle('Velocità')

%velocità in y
vy=cumtrapz(acc(:,2),t);
newVy=cumtrapz(newAcc(:,2),t);

%grafico velocità in y
figure(Name="Asse Y")
subplot(1,2,1)
plot(t,acc(:,2))
grid on
subtitle('Accelerazione')
subplot(1,2,2)
plot(t,vy(:))
grid on
subtitle('Velocità')

%velocità in z
vz=cumtrapz(acc(:,3),t);
newVz=cumtrapz(newAcc(:,3),t);

%grafico velocità in z
figure(Name="Asse Z")
subplot(1,2,1)
plot(t,acc(:,3))
grid on
subtitle('Accelerazione')
subplot(1,2,2)
plot(t,vz(:))
grid on
subtitle('Velocità')

%grafico nuova velocità in x
figure(Name="Asse X")
subplot(1,2,1)
plot(t,newAcc(:,1))
grid on
subtitle('Nuova Accelerazione')
subplot(1,2,2)
plot(t,newVx(:))
grid on
subtitle('Nuova Velocità')

%grafico nuova velocità in y
figure(Name="Asse Y")
subplot(1,2,1)
plot(t,newAcc(:,2))
grid on
subtitle('Nuova Accelerazione')
subplot(1,2,2)
plot(t,newVy(:))
grid on
subtitle('Nuova Velocità')

%grafico nuova velocità in z
figure(Name="Asse Z")
subplot(1,2,1)
plot(t,newAcc(:,3))
grid on
subtitle('Nuova Accelerazione')
subplot(1,2,2)
plot(t,newVz(:))
grid on
subtitle('Nuova Velocità')

% %% Filtro
% %Accelerazione filtrata
% accf=lowpass(newAcc, 40, 1e3);
% %plot
% figure(Name="Accelereazione Filtrata")
% subplot(3,1,1);
% plot(t,accf(:,1),LineWidth=1,Color="b");
% subtitle("X", Color="b");
% subplot(3,1,2);
% plot(t,accf(:,2),LineWidth=1,Color="r");
% subtitle("Y", Color="b");
% subplot(3,1,3);
% plot(t,accf(:,3),LineWidth=1,Color="g");
% subtitle("Z", Color="b");
% 
% %Velocità filtrata
% %velocità in x
% vfx=cumtrapz(t,accf(:,1));
% 
% %grafico velocità filtrata in x
% figure(Name="Asse X")
% subplot(1,2,1)
% plot(t,accf(:,1))
% grid on
% subtitle('Accelerazione filtrata')
% subplot(1,2,2)
% plot(t,vfx(:))
% grid on
% subtitle('Velocità filtrata')
% 
% %velocità filtrata in y
% vfy=cumtrapz(accf(:,2),t);
% 
% %grafico velocità filtrata in y
% figure(Name="Asse Y")
% subplot(1,2,1)
% plot(t,accf(:,2))
% grid on
% subtitle('Accelerazione filtrata')
% subplot(1,2,2)
% plot(t,vfy(:))
% grid on
% subtitle('Velocità filtrata')
% 
% %velocità filtrata in z
% vfz=cumtrapz(accf(:,3),t);
% 
% %grafico velocità filtrata in y
% figure(Name="Asse Z")
% subplot(1,2,1)
% plot(t,accf(:,3))
% grid on
% subtitle('Accelerazione filtrata')
% subplot(1,2,2)
% plot(t,vfz(:))
% grid on
% subtitle('Velocità filtrata')
% 
