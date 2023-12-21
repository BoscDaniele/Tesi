function [gRot,gMedio] = GRot(inizio, fine, dbg)

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
end