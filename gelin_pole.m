clear
clc
%% input
data=[2015;10;24;20;00;00];%% year, month, day, hour, minute,second
out0=[25498.66;0.00456;136.5;63.9;99.7;0];% a(km),e,RA(degrees),incl(degrees),w(degrees),TA(degrees)
step=600; %%seconds
day=30; %% day
%%
S0=time(data);
JD=Jul(data);MJD=JD-2400000.5;
XX0=kparameter(out0);
[X,Y,Z,Vx,Vy,Vz]=trans1(XX0,S0);%%
wz=7.292115e-5;
a=0;b=0;
VV=[Vx;Vy;Vz];
rr=[X;Y;Z];
OM=[0,wz,0;-wz,0,0;0,0,0];
mu=398600.4328969391;
day=day;
i=0;
itt=0:step/86400:day;
%%
for it=0:step/86400:day
    JD1=JD+it;
    i=i+step;
S=(S0+(360.98564724)*i/3600/24)*pi/180;
keps=S*180/pi;
dxy=polarMotion(JD1-2400000.5);
dx=dxy(1);dy=dxy(2);

OM=wz*[0,1,dy;-1,0,dx;-dy,-dx,0];

V=[Vx;Vy;Vz];
h=step;
%-----------------------------
k11=[Vx;Vy;Vz];
r=[X;Y;Z];
k21=-mu/(norm(r))^3*r-OM^2*r+2*OM*k11;
k12=k11+h/2*k21;
r1=r+h/2*(k11);
k22=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k12;
k13=k11+h/2*k22;
r1=r+h/2*(k12);
k23=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k13;
k14=k11+h*k23;
r1=r+h*(k13);
k24=-mu/(norm(r1))^3*r1-OM^2*r1+2*OM*k14;
r=r+h/6*(k11+2*(k12)+2*(k13)+k14);
X=r(1);
Y=r(2);
Z=r(3);
V=V+h/6*(k21+2*k22+2*k23+k24);
Vx=V(1);
Vy=V(2); 
Vz=V(3);

%-----------------------------
a=a+1;
VVX(a)=Vx;
VVY(a)=Vy;
VVZ(a)=Vz;
Rx(a)=X;
Ry(a)=Y;
Rz(a)=Z;
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
legend('?X','?Y','?Z')
set(gca,'linewidth',1,'fontsize',16);

figure(2)
subplot(2,3,1)
plot(itt,kepler_a,'k','LineWidth',1)
grid on
hold on
ylabel('км');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('большая полуось')

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
ylabel('градус');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('наклонение')

subplot(2,3,4)
plot(itt,kepler_RAAN,'r','LineWidth',1)
grid on
hold on
ylabel('градус');
xlabel('день');
set(gca,'linewidth',1,'fontsize',12);
title('долгота восходящего узла')

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