function [C] = protinani_z_uhlu(Ay,Ax,By,Bx,omegaB,omegaA,poz,uh,slop)
%Výpočet protínání úhlů v pravoúhlé soustavě
%vstup:
%    Ay-souřadnice bodu A Y v metrech
%    Ax-souřadnice bodu A X v metrech
%    By-souřadnice bodu B Y v metrech
%    Bx-souřadnice bodu B X v metrech
%omegaA-úhel u vrcholu A
%omegaB-úhel u vrcholu B
%   poz-pozice bodu C v závislosti spojnice bodu A-B
%    uh-úhel, typ úhlu co je používán 1=radiány (základně nastavené)
%       2=stupně (může být zadáno v [° ' "]) a 3=Gony
%  slop-seřazení jak má být seřazen výsledek 1=[Y,X] (základně nastaven),
%       2=[X,Y] jsou v metrech
%výstup
%     C-Výsledné souřadnice bodu C
%%
    if nargin==7
        uh=1;
    end
    if nargin==8
        slop=1;
    end

    st2rad=pi/180;
    G2rad=pi/200;
    
    if uh==2
        if length(omegaB)~=1
            [omegaB]=dms2deg(omegaB);
        end
        omegaB=omegaB*st2rad;
        if length(omegaA)~=1
            [omegaA]=dms2deg(omegaA);
        end
        omegaA=omegaA*st2rad;
    elseif uh==3
        omegaA=omegaA*G2rad;
        omegaB=omegaB*G2rad;   
    end
    
    smAB=atan2(By-Ay,Bx-Ax);
    if smAB<0
        smAB=smAB+2*pi;
    end
    smBA=smAB-pi;
    if smBA<0
        smBA=smBA+2*pi;
    end
    s=sqrt((Ax-Bx)^2+(Ay-By)^2);
    
    wC=pi-(omegaA+omegaB);
    
    Sac=(s/sin(wC))*sin(omegaB);
    Sbc=(s/sin(wC))*sin(omegaA);
    
    if poz==1
        smAC=smAB+omegaA;
        smBC=smBA-omegaB;
    else
        smAC=smAB-omegaA;
        smBC=smBA+omegaB;
    end
    
    CA=[Ay+Sac*sin(smAC),Ax+Sac*cos(smAC)];
    CB=[By+Sbc*sin(smBC),Bx+Sbc*cos(smBC)];
    roz=(CA-CB);
    C=[mean([CA(1),CB(1)]),mean([CA(2),CB(2)])];
    
    if slop==2
        C=[C(2),C(1)];
    end
    
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