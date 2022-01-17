function [X,Y,Z,Vx,Vy,Vz]=trans1(XX0,S0)
XX=[XX0(1);XX0(2);XX0(3)];
VV=[XX0(4);XX0(5);XX0(6)];
S=S0;
S=S*pi/180;
wz=7.292115e-5;
OM=[0,wz,0;-wz,0,0;0,0,0];
A=[cos(S),sin(S),0;-sin(S),cos(S),0;0,0,1];
XX=A*XX;
VV=A*VV+OM*XX;
X=XX(1);Y=XX(2);Z=XX(3);
Vx=VV(1);Vy=VV(2);Vz=VV(3);