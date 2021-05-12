function StrDis=EigStr_general(Slowness, Polarization, Cij0, freq_angular)
%% Documentation
% Function to compute the stress-displacement vector
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
global p;
Tau=zeros(3,6);
Cij=squeeze(Cij0);

for i=1:6
    k=freq_angular*[p;0;Slowness(i)];
    Strain=1i*[Polarization(1,i)*k(1);...
            0;...
            Polarization(3,i)*k(3);...
            Polarization(2,i)*k(3);...
            Polarization(1,i)*k(3)+Polarization(3,i)*k(1);...
            Polarization(2,i)*k(1)];
    Tau0=Cij*Strain;
    Tau(:,i)=[Tau0(5); Tau0(4); Tau0(3)];
end
    StrDis=[Polarization;Tau];



