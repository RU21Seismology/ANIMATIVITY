function [Slowness_hs, Polarization_hs]=EigHs(para_hs)
%% Documentation
% Function to slowness and polarizations of the half space 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
global p velocity_n
a=para_hs(1)/velocity_n;
b=para_hs(2)/velocity_n;

ra=sqrt(1/(a*a)-p*p);% P wave
rb=sqrt(1/(b*b)-p*p);% SV wave and SH wave

%Allocate memory for each eigen value at each time point
Slowness_hs=[ra; rb; rb; -ra; -rb; -rb];
Polarization_hs=[a*p, b*rb, 0, a*p, b*rb, 0;...
                 0,   0,    1,  0,   0,   -1;...
                 a*ra, -b*p,0, -a*ra, b*p, 0];

        
 