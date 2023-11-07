function [rho,alpha] = CalcolaAngolo(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rho=sqrt(x^2+y^2);

alpha = acos(x/rho);
end