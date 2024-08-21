clc; clear; format long G
global G2R;
G2R=pi/180;

% A = [ zemep. sirka,  zemep. delka, azimut] (s.s.,v.d.)
A = [ 49 56  7.691;   14 45 32.263];
fi=A(1,:); lambda=A(2,:);
[fi1]=dms2deg(fi); fi=90-fi1;
[lambda]=dms2deg(lambda);
AZ= [252 40 34.426];
Pr= [0 0 0]; [Pr]=dms2deg(Pr); Pr=90-Pr;
[Az]=dms2deg(AZ); Az=360-Az;
var=["φ1";"λ1";"Az1";"φ2";"λ2";"Az2"];

%% Příklad 1:
%% A) Prusecik s rovnikem
Rov=[Pr,fi,Az];
[Rov_vys,Rov1] = spher_trig(Rov,'SSu');
Az1=Rov1(5)+180;Az2=360+180-(Rov1(5)+180);
Lam1=lambda-Rov1(6);
Lam2=Lam1+180;
PRUS_rov=[90-Pr;Lam1;Az1;90-Pr;Lam2;Az2];[PRUS_rov]=deg2dms(PRUS_rov);PRUS_rov(:,3)=round(PRUS_rov(:,3),1);
Rov_vys=table(var,PRUS_rov,'VariableNames',["Variable","Value"]);
%% B) Prusecik s 0 Polednikem
Pol=[lambda,fi,Az];
[Pol_vys,Pol1] = spher_trig(Pol,'uSu');
Prus_pol=[90-Pol1(1);0;Pol1(5)+180];[Prus_pol]=deg2dms(Prus_pol);Prus_pol(:,3)=round(Prus_pol(:,3),1);
Pol_vys=table(var(1:3),Prus_pol,'VariableNames',["Variable","Value"]);
%% C) Nejiznejsi a Nejsevernejsi
Jiz=[fi,Az,90];
[JS_vys,Jiz1] = spher_trig(Jiz,'Suu');
fiJ=90-Jiz1(3); LamJ=lambda-Jiz1(5);Azj=270;
fiS=-fiJ;LamS=LamJ+180;Azs=270;
J_S=[fiJ;LamJ;Azj;fiS;LamS;Azs];
var1=["φJ";"λJ";"AzJ";"φS";"λS";"AzS"];
[J_S]=deg2dms(J_S);J_S(:,3)=round(J_S(:,3),1);
JS_vys=table(var1,J_S,'VariableNames',["Variable","Value"]);

%% Příklad číslo 2:
var2=["φP1";"λP1";"φP2";"λP2"];
% C = [ zemep. sirka,  zemep. delka] (s.s.,v.d.)
C = [ 49 36 10.329;   16  6  1.935];
% D = [ zemep. sirka,  zemep. delka] (s.s.,v.d.)
D = [ 50 40 59.217;   15 47  8.460];

S=J_S(4:5,:);
[A_xyz] = sph_2_cart(A(1,:),A(2,:),6378000);
[B_xyz] = sph_2_cart(S(1,:),S(2,:),6378000);
[C_xyz] = sph_2_cart(C(1,:),C(2,:),6378000);
[D_xyz] = sph_2_cart(D(1,:),D(2,:),6378000);

n1=cross(A_xyz,B_xyz);
n2=cross(C_xyz,D_xyz);
p=cross(n1,n2);
p=p./norm(p);
p1=-p;
sphr= cart_2_sphr(p);
prusecik1=deg2dms(sphr(1:2)');
sphr1 = cart_2_sphr(p1);
prusecik2=deg2dms(sphr1(1:2)');
prus=[prusecik1;prusecik2];prus(:,3)=round(prus(:,3),1);
Vys_prusecik=table(var2,prus,'VariableNames',["Variable","Value"]);
%% Příklad 3:
Azi_A=dms2deg(AZ);
C1=dms2deg(C);
D1=dms2deg(D);
vstup=[90-C1(1),C1(2)-D1(2),90-D1(1)];
[Tab_Azi_C,Azi_C] = spher_trig(vstup,'SuS');
Azi_C=360-Azi_C(5);
fi=fi1;
for n=1:7
    k=n-4;
    U(n)=2*(atand(exp((tand(Azi_A)*log(tand(fi/2+45))-tand(Azi_C)*log(tand(C1(1)/2+45))+C1(2)*G2R-lambda*G2R+2*k*pi)/(tand(Azi_A)-tand(Azi_C))))-45);
    V(n)=(tand(Azi_A)*(log(tand(U(n)/2+45))-log(tand(fi/2+45)))+lambda*G2R)/G2R;
    Vt(n)=(tand(Azi_C)*(log(tand(U(n)/2+45))-log(tand(C1(1)/2+45)))+C1(2)*G2R)/G2R;
    while V(n)<-180 || V(n)>180
        if V(n)<-180
            V(n)=V(n)+360;
        elseif V(n)>180
            V(n)=V(n)-360;
        end
    end
    while Vt(n)<-180 || Vt(n)>180
        if Vt(n)<-180
            Vt(n)=Vt(n)+360;
        elseif Vt(n)>180
            Vt(n)=Vt(n)-360;
        end
    end
    pom=[90-C1(1),abs(C1(2)-V(n)),90-U(n)];
    [Tab_del_C,del] = spher_trig(pom,'SuS');
    delka(n)=del(3);
end
delka=delka*G2R*6378000;delka=delka';
ind=find(delka==min(delka));
fi_prus_lox=deg2dms(U(ind));fi_prus_lox(:,3)=round(fi_prus_lox(:,3),1);
lambda_prus_lox=deg2dms(V(ind));lambda_prus_lox(:,3)=round(lambda_prus_lox(:,3),1);
Delka=round(delka(ind),1);
Lox_prus=table(fi_prus_lox,lambda_prus_lox,Delka,'VariableNames',["φ","λ","delka [m]"]);
%% Porovnání průsečíků
[Tab_roz_P,roz_P] = spher_trig([90-U(ind),abs(sphr(2)-V(ind)),90-sphr(1)],'SuS');
roz=roz_P(3)*G2R*6378000;