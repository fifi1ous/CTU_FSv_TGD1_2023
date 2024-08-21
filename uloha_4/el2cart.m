function [XYZ,X,Y,Z] = el2cart(Body,a,b,L,H)
%%
%INPUT         B - Elipsoidical coordinates of just B, or can be matrix of [B,L,H] input in [° ' "], or in [°], H is in [m] 
%              a - Semi-major axis of ellipsoid in [m]
%              b - Semi-minor axis of ellispoid in [m]
%              L - Elipsoidical coordinates of just L in format [° ' "]
%              H - Elipsoidical height in [m]
%OUTPUT      XYZ - Geocentrical coordinates  X, Y, Z [m]
%              X - Geocentrical coordinates  X [m]
%              Y - Geocentrical coordinates  Y [m]
%              Z - Geocentrical coordinates  Z [m]
%%
if nargin==5
    Body=[Body,L,H];
end
e=(a^2-b^2)/(a^2);

if size(Body,2) == 7
    B=dms2deg(Body(:,1:3));
    L=dms2deg(Body(:,4:6));
    H=Body(:,end);
else
   B=Body(:,1);
   L=Body(:,2);
   H=Body(:,3); 
end


N=a./sqrt(1-e*sind(B).^2);

X=(N+H).*cosd(B).*cosd(L);
Y=(N+H).*cosd(B).*sind(L);
Z=(N-e*N+H).*sind(B);
XYZ=[X,Y,Z];

function [v]=deg2dms(st)
%převod stupně v desetinné podobě ne matici
%vstup:
%   st-matice stupňů v desetinné podobě
%výstup:
%   v-matice [s m v]
r=length(st);
v=ones(r,3);
for i=1:r
    s=abs(fix(st(i,1)));
    m1=((abs(st(i,1))-abs(fix(st(i,1))))*60);
    m=fix(m1);
    vt=(m1-fix(m1))*60;
    v(i,1)=s;
    v(i,2)=m;
    v(i,3)=vt;
    if st(i,1)<0
        if s~=0
            v(i,1)=s*(-1);
        elseif s==0&&m~=0
            v(i,2)=m*(-1);
        else s==0&&m==0
            v(i,3)=vt*(-1);
        end
    end
end
end
function [st]=dms2deg(v)
%převod stupňů, minut, vteřin na desetiné číslo
%vstup:
%   v-matice [s m v]
%vystup:
%   st-desetiné číslo
[r]=size(v,1);
st=ones(r,1);
for i=1:r
        if any(v(i,:)<0)
            Q=(-1);
        else
            Q=(1);
        end
        st(i,1)=(abs(v(i,1))+abs(v(i,2)/60)+abs(v(i,3)/3600))*Q;
end
end
end