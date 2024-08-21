function [h_vyr,h0,m0,mh,mH,mX,residua,SS_tr] = hel_7_tr_dif(XYZ_BES,XYZ_WGS)
%%
%INPUT   XYZ_BES - kartesian coordinates on ellipsoid XYZ in [m] (in to
%                  which you want to transform)
%        XYZ_WGS - kartesian coordinates on ellipsoid XYZ in [m] (from
%                  which you want to transform)
%OUTPUT    h_vyr - transformation key
%            h_0 - transformation key reducet do the the center
%             m0 - aposterial standart deviation
%             mh - standart deviation of transformation key for reduced
%                  coordinates
%             mH - standart deviation of transformation key
%             mX - standart deviation of transformed coordinates
%        residua - residue of transformed coordinates
%          SS_tr - transformed coordinates used in key
%%
X_TW=mean(XYZ_WGS(:,1));
Y_TW=mean(XYZ_WGS(:,2));
Z_TW=mean(XYZ_WGS(:,3));

XYZ_W_red=[X_TW-XYZ_WGS(:,1),Y_TW-XYZ_WGS(:,2),Z_TW-XYZ_WGS(:,3)];

X_TB=mean(XYZ_BES(:,1));
Y_TB=mean(XYZ_BES(:,2));
Z_TB=mean(XYZ_BES(:,3));

XYZ_B_red=[X_TB-XYZ_BES(:,1),Y_TB-XYZ_BES(:,2),Z_TB-XYZ_BES(:,3)];

h0=[0;0;0;1;0;0;0];
xs2=[]; xs1=[];
for n=1:size(XYZ_W_red,1)
    xs2=[xs2;XYZ_B_red(n,1);XYZ_B_red(n,2);XYZ_B_red(n,3)];
    xs1=[xs1;XYZ_W_red(n,1);XYZ_W_red(n,2);XYZ_W_red(n,3)];
end

o=0;
v2=0;v1=1;

l0=trsnsformace(XYZ_W_red,h0);
while any(abs(v2-v1)>0.0001)
    [A] = matice_A(XYZ_W_red,h0);
    lc=l0-xs2;
    dh=-(A'*A)^(-1)*A'*lc;
    h0=h0+dh;
    
    v1=A*dh+lc;
    l0=trsnsformace(XYZ_W_red,h0);
    v2=l0-xs2;
    o=o+1;
end


m0=sqrt((v1'*v1)/(length(xs1)-7));
Q=(A'*A)^(-1);
mh=m0*sqrt(diag(Q));
Qx=A*(A'*A)^(-1)*A';
mx=m0*sqrt(diag(Qx));

Translace=[X_TB;Y_TB;Z_TB]-h0(1:3)-h0(4)*[1,h0(7),-h0(6);-h0(7),1,h0(5);h0(6),-h0(5),1]*[X_TW;Y_TW;Z_TW];
h_vyr=[Translace;h0(4:end)];
SS_tr1=trsnsformace(XYZ_WGS,h_vyr);

mX=[];
residua=[];
SS_tr=[];

A_pom=matice_A([X_TW,Y_TW,Z_TW],h0);
A_pom=A_pom(:,4:end).^2;
Pr_T=sqrt(eye(3)*mh(1:3).^2+A_pom*mh(4:end).^2);
mH=[Pr_T;mh(4:end)];

for n=1:3:length(xs2)
    mX=[mX;mx(n:n+2)'];
    residua=[residua;v2(n:n+2)'];
    SS_tr=[SS_tr;SS_tr1(n:n+2)'];
end

function [A] = matice_A(body,T)
TA=eye(3,3);
alpha=T(5);
beta=T(6);
gama=T(7);
q=T(4);

A=[];
for n=1:size(body,1)
    q1=body(n,1)+gama*body(n,2)-beta*body(n,3);
    q2=-gama*body(n,1)+body(n,2)+alpha*body(n,3);
    q3=beta*body(n,1)-alpha*body(n,2)+body(n,3);
    R=[0,               -q*body(n,3),   q*body(n,2);
       q*body(n,3),     0,              -q*body(n,1);
       -q*body(n,2),    q*body(n,1),    0];
    A1=[TA,[q1;q2;q3],R];
    A=[A;A1];
end
end

function [lvyr] = trsnsformace(body,T)
    TA=[T(1);T(2);T(3)];
    q=T(4);
    R=[1,T(7),-T(6);-T(7),1,T(5);T(6),-T(5),1];
    
    lvyr=[];
    for n=1:length(body)
        lvyr=[lvyr;TA+q*R*[body(n,1);body(n,2);body(n,3)]];
    end
end

end