function S0=time(data)
year=data(1);
month=data(2);
day=data(3);
hour=data(4);
minute=data(5);
second=data(6);
L=0;

j0=367*year-fix(7*(year+fix((month+9)/12))/4)+fix(275*month/9)+day+1721013.5;
ut=hour+minute/60+second/3600;
T0=(j0-2451545)/36525;
theta0=100.4606184+36000.77004*T0+0.000387933*T0^2-2.583e-8*T0^3;
if theta0/360>1
    theta0=theta0-fix(theta0/360)*360;
end
thetaG=theta0+360.98564724*ut/24;
theta=thetaG+L;


if theta/360>1
    theta=theta-fix(theta/360)*360;
end

S0=theta;
