function [xyz] = sph_2_cart(fi,lambda,ro)
%%
%INPUT      fi     - elevation, mandatory. Either single value, or all values
%                    of fi lambda . Matrix of Degrees, minutes and
%                    seconds or vector of decimal numbers
%           lambda - azimuth. As a vector or matrix 
%           ro     - radius of spehere, if not used it will calculet as
%                    sphere with radius 1
%OUTPUT     Cartesian - xyz, vector of cartesian coordinates 
%%
if nargin==1
    lambda=fi(2,:);
    fi=fi(1,:);
    ro=1;
end
if nargin ==2
    ro=1;
end
if size (fi,2)>1
    [fi]=dms2deg(fi);
end
if size(lambda,2)>1
    [lambda]=dms2deg(lambda);
end
xyz(1)=ro*cosd(fi)*cosd(lambda);
xyz(2)=ro*cosd(fi)*sind(lambda);
xyz(3)=ro*sind(fi);
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