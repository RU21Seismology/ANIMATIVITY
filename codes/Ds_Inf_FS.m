function Rud0=Ds_Inf_FS(E1, w, z1, z_fs, D)
%% Documentation
% Function to compute the displacement-stress matrix at the free surface 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
E1_21=E1(4:6,1:3);
E1_22=E1(4:6,4:6);

% Build up the reflection and transmission coefficients
% interface 1
delta2_u_z=[exp(-1i*w*D(4)*(z_fs-z1)) exp(-1i*w*D(5)*(z_fs-z1)) exp(-1i*w*D(6)*(z_fs-z1))];

deltaD_RT=diag(delta2_u_z);
Rud0=-E1_21\E1_22*deltaD_RT;




