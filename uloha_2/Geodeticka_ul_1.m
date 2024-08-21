function [FiB,LamdaB,AlphaB,n,h] = Geodeticka_ul_1(FiA,LamA,AziA,sAB,mez,a,e2)
%%
% INPUT       FiA - Latitude of the point A
%            LamA - Longtitude of the point A
%            AziA - Azimuth in the point A
%             sAB - Length in meters of the trendline
%             mez - Precision of the calculation
%               a - Main half-axis of the Ellipsoid
%              e2 - Second excentricity of the Ellipsoid

% OUTPUT      FiB - Latitude of the point B
%            LamB - Longtitude of the point B
%            AziA - Azimuth in the point B
%               n - number of steps
%               h - step in meters

%%
    [mez]=dms2deg(mez);
    mez=mez/180*pi;
    n=1;
    FiB=0;fi=10000;LambdaB=0;lambda=10000;AlphaB=0;alpha=10000;
    AB=sAB;
    while(abs(FiB-fi)>mez||abs(LamdaB-lambda)>mez||abs(AlphaB-alpha)>mez)
        FiB=fi;
        LamdaB=lambda;
        AlphaB=alpha;
    
        fi=FiA/180*pi;
        lambda=LamA/180*pi;
        alpha=AziA/180*pi;
    
        h=AB/n;
        for m=1:n    
            [kif,kalpha,klambda]=koeficienty(h,fi,alpha,a,e2);
            fi=fi+(1/6)*h*(kif(1)+2*kif(2)+2*kif(3)+kif(4));
            lambda=lambda+(1/6)*h*(klambda(1)+2*klambda(2)+2*klambda(3)+klambda(4));
            alpha=alpha+(1/6)*h*(kalpha(1)+2*kalpha(2)+2*kalpha(3)+kalpha(4));
        end
    
        if n==1
            n=n+1;
        else
            n=n+2;
        end
    end
    [FiB]=deg2dms(fi/pi*180);
    [LamdaB]=deg2dms(lambda/pi*180);
    [AlphaB]=deg2dms(alpha/pi*180);
    
    
    function [fiB]=f1(fi,alpha,a,e2)
        M=a*(1-e2)/((sqrt(1-e2*sin(fi)^2))^3);
        fiB=cos(alpha)/M;
    end
    function [LamB]=f2(fi,alpha,a,e2)
        N=a/(sqrt(1-e2*sin(fi)^2));
        LamB=sin(alpha)/(N*cos(fi));
    end
    function [AlphaB]=f3(fi,alpha,a,e2)
        N=a/(sqrt(1-e2*sin(fi)^2));
        AlphaB=(sin(alpha)*tan(fi))/N;
    end
    function [kif,kalpha,klambda]=koeficienty(h,fi,alpha,a,e2)
        kif(1)=f1(fi,alpha,a,e2);                                   kalpha(1)=f3(fi,alpha,a,e2);
        kif(2)=f1(fi+h/2*kif(1),alpha+h/2*kalpha(1),a,e2);          kalpha(2)=f3(fi+h/2*kif(1),alpha+h/2*kalpha(1),a,e2);
        kif(3)=f1(fi+h/2*kif(2),alpha+h/2*kalpha(2),a,e2);          kalpha(3)=f3(fi+h/2*kif(2),alpha+h/2*kalpha(2),a,e2);
        kif(4)=f1(fi+h*kif(3),alpha+h*kalpha(3),a,e2);              kalpha(4)=f3(fi+h*kif(3),alpha+h*kalpha(3),a,e2);
    
        klambda(1)=f2(fi,alpha,a,e2);                          
        klambda(2)=f2(fi+h/2*kif(1),alpha+h/2*kalpha(1),a,e2); 
        klambda(3)=f2(fi+h/2*kif(2),alpha+h/2*kalpha(2),a,e2); 
        klambda(4)=f2(fi+h*kif(3),alpha+h*kalpha(3),a,e2); 
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