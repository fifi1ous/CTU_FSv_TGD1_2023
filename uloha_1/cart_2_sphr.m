function [sphr,sphrtab] = cart_2_sphr(x,y,z)
%%
%INPUT         x - cartesian coordinates of x-axis, mandatory value in this variable can be stored X,Y,Z 
%              y - cartesian coordinates of y-axis
%              z - cartesian coordinates of z-axis
%OUTPUT     sphr - vector of φ λ ρ values in decimal degrees
%        sphrtab - table of φ λ ρ values in decimal degrees
%%
if nargin==1
    y=x(2);
    z=x(3);
    x=x(1);
end
lambda=atan2d(y,x);
fi=atand(z/sqrt(x^2+y^2));
ro=sqrt(x^2+y^2+z^2);
sphr=[fi,lambda,ro];
if nargout >1
    sphrtab=table(fi,lambda,ro,'VariableNames',["φ","λ","ρ"]);
end
end