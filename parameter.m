function out2=parameter(out1)
mu=398600.4328969391;
eps=1.e-10;
R0x=out1(1);
R0y=out1(2);
R0z=out1(3);
V0x=out1(4);
V0y=out1(5);
V0z=out1(6);

R=[R0x R0y R0z];
V=[V0x V0y V0z];

r=norm(R);
v=norm(V);
vr=dot(R,V)/r;
H=cross(R,V);
h=norm(H);
incl=acos(H(3)/h);
N=cross([0 0 1],H);
n=norm(N);
if n~=0
RA=acos(N(1)/n);
if N(2)<0
RA=2*pi-RA;
end
else
RA=0;
end
E=1/mu*((v^2-mu/r)*R-r*vr*V);
e=norm(E);
if n~=0
	if e>eps
		w=acos(dot(N,E)/n/e);
		if E(3)<0
	w=2*pi-w;
	end
else
	w=0;
end	
else
	w=0;
end

if e>eps
	TA=acos(dot(E,R)/e/r);
	if vr<0
		TA=2*pi-TA;
	end
else
	cp=corss(N,R);
	if cp(3)>=0
		TA=acos(dot(N,R)/n/r);
else
		TA=2*pi-acos(dot(N,R)/n/r);
	end
end

a=h^2/mu/(1-e^2);
RA=RA*180/pi;
incl=incl*180/pi;
w=w*180/pi;
TA=TA*180/pi;
coe=[ h e RA incl w TA a];
out2=coe;
