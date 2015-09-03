%This subroutine uses equations from Reference 1 to determine the pressure
%(in psi) at a specific depth and latitude (lat).

function [p_psi]=find_pres(depth,lat)
%depth - units of meters with depth positive downward
%lat - units of degrees

syms p real     %looking for only real solutions for pressure (p)

%entering in the constants for Eqn 25 from Reference 1 
%and listed in 01 SwRI Rpt No 18-11975-1-Rev-0 - dtd 2-20-06.pdf
C1=9.72659;
C2=-2.2512E-5;
C3=2.279E-10;
C4=-1.82E-15;
gamma=2.184e-6; %mean vertical gradient of gravity, m/s^2/dbar

%assumption from 01 SwRI Rpt No 18-11975-1-Rev-0 - dtd 2-20-06.pdf
deltaD=0;

%finding angle in radians from latitude input and value for gravitational
%constant glat
angle=lat/180*pi;
glat=9.780318*(1+(5.2788e-3)*(sin(angle))^2+(2.36e-5)*(sin(angle))^4);

%solving for pressure (p) at specified depth and latitude
soln=solve((C1*p+C2*p^2+C3*p^3+C4*p^4)/(glat+0.5*gamma*p)+deltaD/9.8-depth,p);

%from different solutions for pressure given fourth-order equation, find
%value of pressure that is greater than 0
for i=1:1:size(soln,1)
    if double(soln(i))>0
        p_dbar=double(soln(i));
    end
end

%convert dbar to psi, where 1 psi = 0.689475728 dbar
p_psi=p_dbar/0.689475728;