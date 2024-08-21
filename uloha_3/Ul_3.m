% TG1 ul3 - zadani c.17
% delkove veliciny: metry
% uhlove veliciny: DD - stupne, MM - minuty, SS.SSSS - vteriny
% ------------------------------------------------------------

% bod A
% Afl = [ zemep. sirka , zemep. delka ] (elipsoidicke souradnice fi,lambda)
  Afl = [ 50  9 17.3836  14 18 12.9053];
% Auv = [ zemep. sirka , zemep. delka ] (sfericke souradnice U,V)
  Auv = [ 50  6 50.0298  31 59 21.6731];
% Asd = [kartog. sirka , kartog. delka] (kartograficke souradnice S,D)
  Asd = [ 78 40 51.0697  36 41  6.7953];
% Are = [   pruvodic   , polarni uhel ] (polarni souradnice ro,epsilon)
  Are = [  1277900.453   35 56 55.5154];
% Axy = [      x       ,      y       ] (kartezske souradnice X,Y)
  Axy = [  1034514.586     750206.064];

% mereny uhel omegaA na bode A
  omegaA = [  48 23 51.0702]; [omegaA]=dms2deg(omegaA);

% bod B
% Bfl = [ zemep. sirka , zemep. delka ] (elipsoidicke souradnice fi,lambda)
  Bfl = [ 49 51 48.8037  15 28  9.2074];
% Buv = [ zemep. sirka , zemep. delka ] (sfericke souradnice U,V)
  Buv = [ 49 49 22.8404  33  9 20.4825];
% Bsd = [kartog. sirka , kartog. delka] (kartograficke souradnice S,D)
  Bsd = [ 78 45 17.3200  32 35 20.7388];
% Bre = [   pruvodic   , polarni uhel ] (polarni souradnice ro,epsilon)
  Bre = [  1269664.891   31 56  5.4903];
% Bxy = [      x       ,      y       ] (kartezske souradnice X,Y)
  Bxy = [  1077501.150     671595.270];

% mereny uhel omegaB na bode B
  omegaB = [  11 32 10.8110]; [omegaB]=dms2deg(omegaB);

% poloha bodu C
  polohaC = -1;
% (polohaC = 1 pro bod C VPRAVO od spojnice z A do B)
% (polohaC = -1 pro bod C VLEVO od spojnice z A do B)
%%
R=6380703.6105;
DAC=0;DBC=0;
dcy=1;dcx=1;
n=1;
[DAB,DBA,delAB,delBA] = smerova_korekce(Axy,Bxy,Are(1,2:4),Bre(1,2:4),Bre(1,1),Are(1,1),Asd(1,1:3),Bsd(1,1:3));
Cvys=[dcx,dcy];
n=1;

while(abs(dcy)>0.001||abs(dcx)>0.001)
    wA(n)=omegaA-DAC+DAB;
    wB(n)=omegaB+DBC-DBA;
    [C] = protinani_z_uhlu(Axy(2),Axy(1),Bxy(2),Bxy(1),wB(n),wA(n),polohaC,2,2);
    [sd, re] = xy2sd(C);
    [DAC,DCA,delAC,delCA] = smerova_korekce(Axy,C,Are(1,2:4),re(2),re(1),Are(1,1),Asd(1,1:3),sd(1));
    [DBC,DCB,delBC,delCB] = smerova_korekce(Bxy,C,Bre(1,2:4),re(2),re(1),Bre(1,1),Bsd(1,1:3),sd(1));
    Cvys=[Cvys;C];
    dcy=Cvys(n,2)-C(2);
    dcx=Cvys(n,1)-C(1);
    n=n+1;
end
naz=["dAB";"dBA";"dAC";"dCA";"dBC";"dCB"];
smerove_korekce=[round(delAB,3);round(delBA,3);round(delAC,3);round(delCA,3);round(delBC,3);round(delCB,3)];
C=[round(C(1),3),round(C(2),3)];
vys_korekce=table(naz,smerove_korekce);

Sab=sqrt((Axy-Bxy)*(Axy-Bxy)');
Sac=sqrt((Axy-C)*(Axy-C)');
Scb=sqrt((Bxy-C)*(Bxy-C)');
S=(Sab+Sac+Scb)/2;
S=sqrt(S*(S-Sab)*(S-Sac)*(S-Scb));
eps11=(180/pi)*3600*(S/(R^2));
eps22=(delAC+delCB+delBA)-(delAB+delBC+delCA);

WA=deg2dms(wA(end));
WB=deg2dms(wB(end));

%% Korekce, Azimutu
smAB=atan2(Bxy(2)-Axy(2),Bxy(1)-Axy(1))/pi*180;
if smAB<0
    smAB=smAB+360;
end
smBA=smAB-180;
if smBA<0
    smBA=smBA+360;
end

Uq=[59,42,42.69689];    [Uq]=dms2deg(Uq);
Vq=[42,31,31.41725];    [Vq]=dms2deg(Vq);

Af=Afl(1:3);[Af]=dms2deg(Af);Al=Afl(4:6);[Al]=dms2deg(Al);
Bf=Bfl(1:3);[Bf]=dms2deg(Bf);Bl=Bfl(4:6);[Bl]=dms2deg(Bl);


GA=asind((sind(dms2deg(Asd(4:end)))*cosd(Uq))/cosd(dms2deg(Auv(1:3))));
CA=dms2deg(Are(2:end))-GA;
AziA=smAB-180-CA-DAB;
AziA=deg2dms(AziA);
convA=deg2dms(CA);

GB=asind((sind(dms2deg(Bsd(4:end)))*cosd(Uq))/cosd(dms2deg(Buv(1:3))));
CB=dms2deg(Bre(2:end))-GB;
AziB=smBA-180-CB-DBA;
AziB=AziB+360;
AziB=deg2dms(AziB);
convB=deg2dms(CB);

%% měřítko
 Q=[59, 42, 42.69689;42,31,31.41725];[Q]=dms2deg(Q);  %souřadnice kartografického pólu Q
 S0=[78,30,0];[S0]=dms2deg(S0);                       %základní kartografická rovnoběžka Š0
 R=6380703.6105;
Sjtsk=sqrt((Axy-Bxy)*(Axy-Bxy)')

 [sd re] = xy2sd(Axy);
 mA=(sind(S0)*re(1))/(R*cosd(sd(1)))

  [sd re] = xy2sd(Bxy);
 mB=(sind(S0)*re(1))/(R*cosd(sd(1)))

 Sxy=[(Axy(1)+Bxy(1))/2,(Axy(2)+Bxy(2))/2];
  [sd re] = xy2sd(Sxy);
 mS=(sind(S0)*re(1))/(R*cosd(sd(1)))

 Sel=Sjtsk*(1/6)*((1/mA)+(4/mS)+(1/mB))
 
