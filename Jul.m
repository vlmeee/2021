function JD=Jul(data)
year=data(1);
month=data(2);
day=data(3);
hour=data(4);
minute=data(5);
second=data(6);

j0=367*year-fix(7*(year+fix((month+9)/12))/4)+fix(275*month/9)+day+1721013.5;
ut=hour+minute/60+second/3600;
jd=j0+ut/24;
JD=jd;



