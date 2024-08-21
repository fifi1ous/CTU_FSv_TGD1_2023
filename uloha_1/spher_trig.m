function [Vys,C1,C,CB] = spher_trig(C,text,p)
%%
% INPUT         C - Matrix or vector of varible
%            text - Type of spherical triangle you want to calculate
%                   'SSS','uuu','SuS','uSu','SSu' or 'Suu', your values
%                   should be sorted counterclockwise
%               p - precision of seconds,volunteraly variable
% OUTPUT      Vys - Table of all values
%              CB - Example of how are the values sorted
%              C1 - Values in vector format, and in decimal format
%               C - Values in matrix format degrees, minutess

%%
if size(C,2)>1
    if size(C,1)==1
        C=C';
    else
        [C]=dms2deg(C);
    end
end
if C(C>=180)
    error("Délka, nebo úhlel je větší než 180°")
end
if size(C,1)<3
    error("Málo údajů nelze dopočítat")
end

CB=["a";"b";"c";"α";"β";"γ"];
if text=='SSS'
    a=C(1);b=C(2);c=C(3);
    [alfa]=kos_str_uhel(a,b,c);
    [beta]=kos_str_uhel(b,a,c);
    [gama]=kos_str_uhel(c,b,a);
elseif text=='uuu'
    alfa=C(1);beta=C(2);gama=C(3);
    [a]=kos_uhel_str(alfa,beta,gama);
    [b]=kos_uhel_str(beta,alfa,gama);
    [c]=kos_uhel_str(gama,beta,alfa);
elseif text=='SuS'
    a=C(1);gama=C(2);b=C(3);
    [c]=kos_str(a,b,gama);
    [alfa]=kos_str_uhel(a,b,c);
    [beta]=kos_str_uhel(b,a,c);
elseif text=='uSu'
    gama=C(1);b=C(2);alfa=C(3);
    [beta]=kos_uhel(gama,alfa,b);
    [a]=kos_uhel_str(alfa,beta,gama);
    [c]=kos_uhel_str(gama,beta,alfa);
elseif text=='SSu'
    a=C(1);b=C(2);alfa=C(3);
    [beta]=sin_uhel(b,a,alfa);
    [beta]=rozpoznani_uhlu(b,a,beta,alfa);    
    if a==90 && b== 90
        if b == 90
            c=90;
        else
            [c]=kos_sin_del(b,alfa);
        end
    else
        [c]=neper_strana(a,b,alfa,beta);
    end
    if alfa ==90 && beta == 90
        if alfa==90
            gama=90;
        else
            [gama]=kos_sin_uh(b,alfa);
        end
    else
        [gama]=neper_uhel(a,b,alfa,beta);
    end
elseif text=='Suu'
    a=C(1);gama=C(2);alfa=C(3);
    [c]=sin_del(gama,a,alfa);
    [c]=rozpoznani_uhlu(gama,alfa,c,a);  
    if a==90 && c== 90
        if c == 90
            b=90;
        else
            [b]=kos_sin_del(c,alfa);
        end
    else
        [b]=neper_strana(a,c,alfa,gama);
    end
    if alfa ==90 && gama == 90
        if alfa==90
            beta=90;
        else
            [beta]=kos_sin_uh(c,alfa);
        end
    else
        [beta]=neper_uhel(a,c,alfa,gama);
    end
end

C1=[a;b;c;alfa;beta;gama];
C=deg2dms(C1);
if nargin==3
    C(:,3)=round(C(:,3),p);
end
Vys=table(CB,C(:,1),C(:,2),C(:,3),'VariableNames',["Variable","DEC °"," '"," """]);

function [v]=deg2dms(st)
%převod stupně v desetinné podobě na stupně vteřiny a minuty
%vstup:
%   st-matice stupňů v desetinné podobě
%výstup:
%   v-matice [s m v]
[r,s]=size(st); q=0;
if r<s
    st=st';
    q=1;
end

r=length(st);
v=ones(r,3);
for i=1:r
    s=abs(fix(st(i,1)));
    m1=((abs(st(i,1))-abs(fix(st(i,1))))*60);
    m=fix(m1);
    vt=(m1-fix(m1))*60;
    if round(vt)==60
        vt=0;
        m=0;
        s=s+1;
    end
    if round(vt,2)==0
        vt=0;
    end
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
if q==1
    v=v';
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
function [a]=kos_str(b,c,alfa)
    a=acosd(cosd(b)*cosd(c)+sind(b)*sind(c)*cosd(alfa));
end
function [alfa]=kos_str_uhel(a,b,c)
    alfa=acosd((cosd(a)-cosd(b)*cosd(c))/(sind(b)*sind(c)));
end
function [alfa]=kos_uhel(beta,gama,a)
    alfa=acosd((-1)*cosd(beta)*cosd(gama)+sind(beta)*sind(gama)*cosd(a));
end
function [a]=kos_uhel_str(alfa,beta,gama)
    a=acosd((cosd(alfa)+cosd(beta)*cosd(gama))/(sind(beta)*sind(gama)));
end
function [a]=sin_del(alfa,b,beta)
    a=asind((sind(alfa)/sind(beta))*sind(b));
end
function [alfa]=sin_uhel(a,b,beta)
    alfa=asind((sind(a)/sind(b))*sind(beta));
end
function [gama]=neper_uhel(a,b,alfa,beta)
    gama=2*acotd((cosd(a/2 + b/2)*tand(alfa/2 + beta/2))/cosd(a/2 - b/2));
end
function [c]=neper_strana(a,b,alfa,beta)
    c=2*atand((cosd(alfa/2 + beta/2)*tand(a/2 + b/2))/cosd(alfa/2 - beta/2));
end
function [alfa]=rozpoznani_uhlu(a,b,alfa,beta)
    if a-b<0
        ab=-1;
    else
        ab=1;
    end
    if alfa-beta<0
        afbt=-1;
    else
        afbt=1;
    end
    if ab+afbt==0
        alfa=180-alfa;
    end    
end
function [c]=kos_sin_del(b,alfa)
    c=atand(cosd(alfa)/cotd(b));
end
function [gama]=kos_sin_uh(a,beta)
    gama=atand(-cosd(a)/cotd(beta));
end
end