clear
clc
%% input
data=[2018;10;24;20;00;00];%% year, month, day, hour, minute,second
out0=[25498.66;0.00456;136.5;63.9;99.7;0];% a(km),e,RA(degrees),incl(degrees),w(degrees),TA(degrees)
step=1000; %%seconds
day=30; %% day
%%
S0=time(data);
JD=Jul(data);MJD=JD-2400000.5;
XX0=kparameter(out0);
[X,Y,Z,Vx,Vy,Vz]=trans1(XX0,S0);%%
wz=7.292115e-5;
a=0;b=0;
VV=[Vx;Vy;Vz];rr=[X;Y;Z];
OM=[0,wz,0;-wz,0,0;0,0,0];
mu=398600.4328969391;
day=day;
i=0;
itt=0:step/86400:day;
%%
for it=0:step/86400:day
    i=i+step;
S=(S0+(360.98564724)*i/3600/24)*pi/180;
keps=S*180/pi;
SS=[cos(S),sin(S),0;-sin(S),cos(S),0;0,0,1];
ttt=MJD+i/86400;
tt=ttt;
tt2=tt+2400000.5;
angles= earthNutation(tt2);
%%
T=(tt-51545)/36525;
A1 = 2.650545+2306.083227*T +0.2988499*T^2+0.01801828*T^3-0.000005971*T^4-0.0000003173*T^5;
A2 = -2.650545 +2306.077181*T +1.0927348*T^2+ 0.01826837*T^3 -0.000028596*T^4- 0.0000002904*T^5;
A3 = 2004.191903*T - 0.4294934*T^2 - 0.04182264*T^3 - 0.000007089*T^4- 0.0000001274*T^5;
e=84381.406- 46.836769*T - 0.0001831*T^2+ 0.00200340*T^3- 0.000000576*T^4- 0.0000000434*T^5;
A1=A1/3600*pi/180;
A2=A2/3600*pi/180;
A3=A3/3600*pi/180;
e=e/3600*pi/180;
f1=19.9/1000/3600*pi/180;
f2=9.1/1000/3600*pi/180;
f3=-22.9/1000/3600*pi/180;
B=[1,0,0;0,cos(f1),sin(f1);0,-sin(f1),cos(f1)]*[cos(f2),0,-sin(f2);0,1,0;sin(f2),0,cos(f2)]*[cos(f3),sin(f3),0;-sin(f3),cos(f3),0;0,0,1];
P=[cos(-A2),sin(-A2),0;-sin(-A2),cos(-A2),0;0,0,1]*[cos(A3),0,-sin(A3);0,1,0;sin(A3),0,cos(A3)]*[cos(-A1),sin(-A1),0;-sin(-A1),cos(-A1),0;0,0,1];
f1=-e-angles(2);
f3=angles(1);
f4=e;
N=[1,0,0;0,cos(f1),sin(f1);0,-sin(f1),cos(f1)]*[cos(-f3),sin(-f3),0;-sin(-f3),cos(-f3),0;0,0,1]*[1,0,0;0,cos(f4),sin(f4);0,-sin(f4),cos(f4)];
%%
tt=ttt+0.0000116*1000;
tt2=tt+2400000.5;
angles= earthNutation(tt2);
T=(tt-51545)/36525;
A1 = 2.650545+2306.083227*T +0.2988499*T^2+0.01801828*T^3-0.000005971*T^4-0.0000003173*T^5;
A2 = -2.650545 +2306.077181*T +1.0927348*T^2+ 0.01826837*T^3 -0.000028596*T^4- 0.0000002904*T^5;
A3 = 2004.191903*T - 0.4294934*T^2 - 0.04182264*T^3 - 0.000007089*T^4- 0.0000001274*T^5;
e=84381.406- 46.836769*T - 0.0001831*T^2+ 0.00200340*T^3- 0.000000576*T^4- 0.0000000434*T^5;
A1=A1/3600*pi/180;
A2=A2/3600*pi/180;
A3=A3/3600*pi/180;
e=e/3600*pi/180;
f1=19.9/1000/3600*pi/180;
f2=9.1/1000/3600*pi/180;
f3=-22.9/1000/3600*pi/180;
B=[1,0,0;0,cos(f1),sin(f1);0,-sin(f1),cos(f1)]*[cos(f2),0,-sin(f2);0,1,0;sin(f2),0,cos(f2)]*[cos(f3),sin(f3),0;-sin(f3),cos(f3),0;0,0,1];
P2=[cos(-A2),sin(-A2),0;-sin(-A2),cos(-A2),0;0,0,1]*[cos(A3),0,-sin(A3);0,1,0;sin(A3),0,cos(A3)]*[cos(-A1),sin(-A1),0;-sin(-A1),cos(-A1),0;0,0,1];
f1=-e-angles(2);
f3=angles(1);
f4=e;
N2=[1,0,0;0,cos(f1),sin(f1);0,-sin(f1),cos(f1)]*[cos(-f3),sin(-f3),0;-sin(-f3),cos(-f3),0;0,0,1]*[1,0,0;0,cos(f4),sin(f4);0,-sin(f4),cos(f4)];
SS1=N2*P2-N*P;
SS0=N*P;
%-----------------------------------------------
SS2=SS*SS0;
%---------------------------------------
V=[Vx;Vy;Vz];
h=step;
SS1=SS1/h;
%%
k11=[Vx;Vy;Vz];
r=[X;Y;Z];
k21=-mu/(norm(r))^3*r-OM^2*r+2*OM*k11+2*SS2*SS1*SS'*(k11-OM*r);
k12=k11+h/2*k21;
r1=r+h/2*(k11);
k22=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k12+2*SS2*SS1*SS'*(k12-OM*r1);
k13=k11+h/2*k22;
r1=r+h/2*(k12);
k23=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k13+2*SS2*SS1*SS'*(k13-OM*r1);
k14=k11+h*k23;
r1=r+h*(k13);
k24=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k14+2*SS2*SS1*SS'*(k14-OM*r1);
r=r+h/6*(k11+2*(k12)+2*(k13)+k14);
X=r(1);Y=r(2);Z=r(3);
V=V+h/6*(k21+2*k22+2*k23+k24);
Vx=V(1);Vy=V(2); Vz=V(3);
%-----------------------------
a=a+1;
VVX(a)=Vx;VVY(a)=Vy;VVZ(a)=Vz;
Rx(a)=X;Ry(a)=Y;Rz(a)=Z;
%%
XXz=[X,Y,Z,Vx,Vy,Vz];
out1=trans2(XXz,keps);
out2=parameter(out1);
kepler_a(a)=out2(7);
kepler_e(a)=out2(2);
kepler_RAAN(a)=out2(3);
kepler_incl(a)=out2(4);
kepler_w(a)=out2(5);
kepler_p(a)=out2(7)*(1-out2(2)^2);
%%

OM=[0,wz,0;-wz,0,0;0,0,0];
kk11=VV;
kk21=-mu/(norm(rr))^3*rr-OM^2*rr+2*OM*kk11;
kk12=kk11+h/2*kk21;
rr1=rr+h/2*(kk11);
kk22=-mu/(norm(rr1))^3*rr1-OM^2*rr1+2*OM*kk12;
kk13=kk11+h/2*kk22;
rr1=rr+h/2*(kk12);
kk23=-mu/(norm(rr1))^3*rr1-OM^2*rr1+2*OM*kk13;
kk14=kk11+h*kk23;
rr1=rr+h*(kk13);
kk24=-mu/(norm(rr1))^3*rr1-OM^2*rr1+2*OM*kk14;
rr=rr+h/6*(kk11+2*(kk12)+2*(kk13)+kk14);
VV=VV+h/6*(kk21+2*kk22+2*kk23+kk24);
Xx=rr(1);
Yx=rr(2);
Zx=rr(3);

Vx2=VV(1);
Vy2=VV(2);
Vz2=VV(3);
rr=[Xx;Yx;Zx];
VV=[Vx2;Vy2;Vz2];
%-----------------------------
b=b+1;
VVX2(b)=Vx2;
VVY2(b)=Vy2;
VVZ2(b)=Vz2;
Rx2(b)=Xx;
Ry2(b)=Yx;
Rz2(b)=Zx;
end
%-----------------------------
for i=1:1:length(Rx)
RRR=[Rx(i);Ry(i);Rz(i)];
VVV=[VVX(i);VVY(i);VVZ(i)];
yy=RRR/norm(RRR);
Va=VVV+cross([0;0;wz],RRR);
zz=cross(Va,yy)/norm(cross(Va,yy));
xx=cross(zz,yy);
XXX1(i)=Rx2(i)-Rx(i);
YYY1(i)=Ry2(i)-Ry(i);
ZZZ1(i)=Rz2(i)-Rz(i);

tx=[xx';yy';zz']*[XXX1(i);YYY1(i);ZZZ1(i)];
XXX1(i)=tx(1);
YYY1(i)=tx(2);
ZZZ1(i)=tx(3);
RR1(i)=sqrt(XXX1(i)^2+YYY1(i)^2+ZZZ1(i)^2)*1000;
end
%%

%%
figure(1)
plot(itt,XXX1,'b','LineWidth',1)
hold on
plot(itt,YYY1,'r','LineWidth',1)
hold on
plot(itt,ZZZ1,'k','LineWidth',1)
grid on
hold on
ylabel('км');
xlabel('день');
% set(gca,'Xtick',[0:1:day])
legend('deltaX','deltaY','deltaZ')
set(gca,'linewidth',1,'fontsize',16);

figure(2)
subplot(2,3,1)
plot(itt,kepler_a,'k','LineWidth',1)
grid on
hold on
ylabel('км');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('больша€ полуось')

subplot(2,3,2)
plot(itt,kepler_e,'r','LineWidth',1)
grid on
hold on
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('эксцентриситет')

subplot(2,3,3)
plot(itt,kepler_incl,'b','LineWidth',1)
grid on
hold on
ylabel('градусы');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('наклонение')

subplot(2,3,4)
plot(itt,kepler_RAAN,'r','LineWidth',1)
grid on
hold on
ylabel('градусы');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('долгота восходщ€его узла')

subplot(2,3,5)
plot(itt,kepler_w,'k','LineWidth',1)
grid on
hold on
ylabel('градус');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('аргумент перицентра')

subplot(2,3,6)
plot(itt,kepler_p,'b','LineWidth',1)
grid on
hold on
ylabel('км');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('фокальный параметр')