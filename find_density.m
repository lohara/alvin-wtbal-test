%This subroutine uses equations from Reference 1 to determine the density
%(in lbs/in^3) given the pressure (in psi) at depth.

function [rho]=find_density(pressure)

%convert psi to dbar, where 1 psi = 0.689475728 dbar and divide by 10 to
%convert to atm
P=pressure*0.689475728/10;
%assuming S and T are  constant; same assumption used in find_pres
S=35;
%T=19;
T=0;

k0=8.50935e-5;	%another value has been added in brackets 3.47718e-5
k1=-6.12293e-6;
k2=5.2787e-8;
B_w=k0+k1*T+k2*T^2;

m0=-9.9348e-7;
m1=2.0816e-8;
m2=9.1697e-10;
B=B_w+(m0+m1*T+m2*T^2)*S;

h0=3.239908;	%another value has been added in brackets -0.1194975
h1=1.43713e-3;
h2=1.16092e-4;
h3=-5.77905e-7;
A_w=h0+h1*T+h2*T^2+h3*T^3;

i0=2.2838e-3;
i1=-1.0981e-5;
i2=-1.6078e-6;
j0=1.91075e-4;
A=A_w+(i0+i1*T+i2*T^2)*S+j0*S^(3/2);

e0=19652.21;	%another value has been added in brackets -1930.06
e1=148.4206;
e2=-2.327105;
e3=1.360477e-2;
e4=-5.155288e-5;
K_w=e0+e1*T+e2*T^2+e3*T^3+e4*T^4;

f0=54.6746;
f1=-0.603459;
f2=1.09987e-2;
f3=-6.167e-5;
g0=7.944e-2;
g1=1.6483e-2;
g2=-5.3009e-4;
K_0=K_w+(f0+f1*T+f2*T^2+f3*T^3)*S+(g0+g1*T+g2*T^2)*S^(3/2);

K=K_0+A*P+B*P^2;

a0=999.842594; %another value has been added in brackets -28.263737
a1=6.793952e-2;
a2=-9.095290e-3;
a3=1.001685e-4;
a4=-1.120083e-6;
a5=6.536332e-9;
rho_w=a0+a1*T+a2*T^2+a3*T^3+a4*T^4+a5*T^5;

b0=8.24493e-1;
b1=-4.0899e-3;
b2=7.6438e-5;
b3=-8.2467e-7;
b4=5.3875e-9;
c0=-5.72466e-3;
c1=1.0227e-4;
c2=-1.6546e-6;
d0=4.8314e-4;
rho_0=rho_w+(b0+b1*T+b2*T^2+b3*T^3+b4*T^4)*S+(c0+c1*T+c2*T^2)*S^(3/2)+d0*S^2;

rho=rho_0/(1-(P/K));

%convert kg/m^3 to lbs/in^3, where 1 kg/m^3 = 3.6127E-5 lbs/in^3
format long
rho=rho*3.6127e-5;