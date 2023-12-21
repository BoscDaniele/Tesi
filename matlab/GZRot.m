function [gzRot,gMedio] = GZRot(path)
% Questa funzione serve a calcolare la matrice di rotazione da applicare ai
% dati raccolti dal sensore al fine di ruotarli per fa coincidere il
% sistema di riferimento del sensore a quello della bicicletta

% Vengono calcolati gli angoli che il vettore gravità ha rispetto alla
% verticale, una volta trovati si ottiene la matrice di rotazione prima in
% X e successivamente in Y

% La rotazione attorno all'asse Z viene ottenuta inclinando la bicicletta e
% misurando la componente lungo l'asse X della gravità

dbg=importdata(path + "BlueCoin_Log_N000.csv").data;

inizio=125;
fine=length(dbg)-125;

g=dbg(inizio:fine,2:4);
vettG=mean(g);

modYZ=sqrt(vettG(2)^2+vettG(3)^2);
thetaYZ=-acos(-vettG(3)/modYZ);
mRotX=[1,0,0;0,cos(thetaYZ),-sin(thetaYZ);0,sin(thetaYZ),cos(thetaYZ)];

newVettG = vettG*mRotX;

modXz=sqrt(newVettG(1)^2+newVettG(3)^2);
thetaXz=-acos(-newVettG(3)/modXz);
mRoty=[cos(thetaXz),0,sin(thetaXz);0,1,0;-sin(thetaXz),0,cos(thetaXz)];

gRot=mRotX*mRoty;

gMedio=-norm(vettG);

% figure
% plot3([0, vettG(1)],[0, vettG(2)],[0, vettG(3)],LineWidth=1,Color="black");
% hold on
% grid
% axis equal
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="g");
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="b");


%% Rotazione attorno all'asse z
dbz=importdata(path + "BlueCoin_Log_N001.csv").data;
inizio=125;
fine=length(dbz)-125;

y=dbz(inizio:fine,2:4);
vettY=mean(y)*gRot;

modY=sqrt(vettY(1)^2+vettY(2)^2);
thetaxy=-acos(vettY(2)/modY);
mRotz=[cos(thetaxy), -sin(thetaxy), 0; sin(thetaxy), cos(thetaxy), 0; 0, 0, 1];

gzRot=gRot*mRotz;

% figure
% plot3([0, vettY(1)],[0, vettY(2)],[0, vettY(3)],LineWidth=1,Color="black");
% hold on
% grid
% axis equal
% plot3([-500,500],[0,0],[0,0],LineWidth=1,Color="r");
% plot3([0,0],[-500,500],[0,0],LineWidth=1,Color="g");
% plot3([0,0],[0,0],[-500,500],LineWidth=1,Color="b");

end