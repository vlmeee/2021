function XX0=kpatameter(out0)
a=out0(1);
e=out0(2);
RA=out0(3);
RA=RA*pi/180;
incl=out0(4);
incl=incl*pi/180;
w=out0(5);
w=w*pi/180;
TA=out0(6);
TA=TA*pi/180;
mu=398600;
h=sqrt(a*mu*(1-e^2));

rp=(h^2/mu)*(1/(1+e*cos(TA)))*(cos(TA)*[1;0;0]+sin(TA)*[0;1;0]);
vp=(mu/h)*(-sin(TA)*[1;0;0]+(e+cos(TA))*[0;1;0]);

R3_W=[cos(RA) sin(RA) 0;-sin(RA) cos(RA) 0; 0 0 1];
R1_i=[1 0 0;0 cos(incl) sin(incl);0 -sin(incl) cos(incl)];
R3_w=[cos(w) sin(w) 0;-sin(w) cos(w) 0; 0 0 1];

Q_pX=R3_W'*R1_i'*R3_w';
r=Q_pX*rp;
v=Q_pX*vp;
r=r';
v=v';

XX0=[r,v];

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
