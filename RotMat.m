function [mRot,mRotY,mRotx] = RotMat(vet)
%calcolo rotazione attorno all'asse Y
modXZ=sqrt(vet(1)^2+vet(3)^2);
thetaXZ=acos(vet(3)/modXZ);

mRotY=[cos(thetaXZ), 0, sin(thetaXZ); 0, 1, 0; -sin(thetaXZ), 0, cos(thetaXZ)].*-1;

%calcolo rotazione attorno al nuovo asse x
vetY=vet*mRotY;

modYz=sqrt(vetY(2)^2+vetY(3)^2);
thetaYZ=acos(vetY(3)/modYz);

mRotx=[1, 0, 0; 0, cos(thetaYZ), -sin(thetaYZ); 0, sin(thetaYZ), cos(thetaYZ)].*-1;

%calcolo matrice di rotazione completa
mRot=mRotY*mRotx;
end