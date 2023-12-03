function [mRot] = RotMat(ang)
ang=ang*-1;

%% rotazione attorno a X
mRotX=[1, 0, 0; 0, cos(ang(1)), -sin(ang(1)); 0, sin(ang(1)), cos(ang(1))];

%% rotazione attorno a Y
mRotY=[cos(ang(2)), 0, sin(ang(2)); 0, 1, 0; -sin(ang(2)), 0, cos(ang(2))];

%% rotazione attorno a Z
mRotZ=[cos(ang(3)), -sin(ang(3)), 0; sin(ang(3)), cos(ang(3)), 0; 0, 0, 1];

%% rotazione totale
mRot=mRotX*mRotY*mRotZ;

end