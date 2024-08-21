clc; clear; format long G
%%
a=6377397.15508;
b=6356078.96290;
%% Zadání
fiA = [ 50 44 30.764];[FiA]=dms2deg(fiA);
lamA=[14 16 36.790  ];[LamA]=dms2deg(lamA);
aziA=[257 13  6.578];[AziA]=dms2deg(aziA);
sAB = 1129562.25;

%% Nastavení
e2=(a^2-b^2)/(a^2);
M=a*(1-e2)/((sqrt(1-e2*sind(FiA)^2))^3);
N=a/(sqrt(1-e2*sind(FiA)^2));
%% Příklad číslo 1
R=sqrt(M*N);
AB=sAB/R/pi*180;
NEZ=[AB,360-AziA,90-FiA];
[PR1_tab,PR1] = spher_trig(NEZ,'SuS',3);

B_PR1=deg2dms([90-PR1(3);LamA-PR1(4);180+PR1(5)]);



%%  Příklad číslo 2
mez=[0,0,0.0001];
[FiBj,LamdaB,AziB] = Geodeticka_ul_1(FiA,LamA,AziA,sAB,mez,a,e2);
B_PR2=[FiBj;LamdaB;AziB];

%%  Príklad číslo 3
% fiB = [ 32 36 6.5378]; 
 [FiB]=dms2deg(FiBj);
% lamB =[ 38 52 56.0979];
 [LamB]=dms2deg(LamdaB);
NEZ=[90-FiA,LamA-LamB,90-FiB];
[PR3_tab,PR3] = spher_trig(NEZ,'SuS');

sab=PR3(3)/180*pi*R;
azia=360-PR3(5);

mez=[0,0,0.00001];
[FiB1,LamdaB1,azib] = Geodeticka_ul_1(FiA,LamA,azia,sab,mez,a,e2);
fib=1000;
lamb=1000;
pom=dms2deg(mez);
% azib=dms2deg(azib);
n=1;

while(abs(fib-FiB)>pom || abs(lamb-LamB)>pom)
[fib]=dms2deg(FiB1);
[lamb]=dms2deg(LamdaB1);
azib=dms2deg(azib);

delfi=(FiB-fib)/180*pi;
dellam=(LamB-lamb)/180*pi;
fistr=((FiB+fib)/2)/180*pi;

N=a/(sqrt(1-e2*sin(fistr)^2));
M=a*(1-e2)/((sqrt(1-e2*sin(fistr)^2))^3);

alphstr=atan2(N*cos(fistr)*dellam,M*delfi);
p=(N*cos(fistr)*dellam)/sin(alphstr);
alphp=alphstr-dellam/2*sin(fistr);

w=azib/180*pi+pi/2-alphp;
dela=(p*cos(w))/sab;
azia=azia+dela/pi*180;

delsab=(((p*cos(w))^2)/(2*sab))+p*sin(w);
sab=sab+delsab;
wv(n)=w;
DSAB(n)=delsab;
DAZI(n)=dela;
IT(n)=n;

[FiB1,LamdaB1,azib] = Geodeticka_ul_1(FiA,LamA,azia,sab,mez,a,e2);
n=n+1;
end

B_PR3=[deg2dms(azia);azib];
SAB_vyp=sab;
B_PR3(:,3)=round(B_PR3(:,3),3);
SAB_vyp=round(SAB_vyp,3);
wv=deg2dms(wv'/180*pi);
DAZI=deg2dms(DAZI'/180*pi);
wv(:,3)=round(wv(:,3),3);
DAZI(:,3)=round(DAZI(:,3),3);

Pomocne_vysledky=[IT',wv,round(DSAB',3),DAZI]