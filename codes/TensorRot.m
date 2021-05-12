function Cij=TensorRot(Cij0,rotdeg)
%% Documentation
% Function to calculate the tensor after rotation 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
% alpha, beta, gamma are the angles that rotates around axes x1, x2, x3
% correspondingly that the three roatation matrices can be written as
% below:

% angle of tilt (°)
a=rotdeg(1)/180*pi;
% angle of dip (°)
b=rotdeg(2)/180*pi;
% angle of back azimuth (°)
c=rotdeg(3)/180*pi;

% rotate around x, a is the tilt angle from vertical 
% positive if CW when look from the above
T1=[cos(a), sin(a), 0;...
    -sin(a), cos(a), 0;...
    0, 0, 1];
% rotate around y, b is the dipping angle from horizontal 
% positive if CW when look from the right side
T2=[1, 0, 0;...
    0, cos(b), sin(b);...
    0, -sin(b), cos(b)];
% rotate around z, c is the azimuth from the North
% positive if CW when look from the above
T3=[cos(c), sin(c), 0;...
    -sin(c), cos(c), 0;...
    0, 0, 1];

% calculate the overall rotation effect
Trot=T3*T2*T1;

% partition the overall rotation into four different sections
T11=Trot.*Trot;
T12=2.*[Trot(1,2)*Trot(1,3), Trot(1,3)*Trot(1,1), Trot(1,1)*Trot(1,2);...
        Trot(2,2)*Trot(2,3), Trot(2,3)*Trot(2,1), Trot(2,1)*Trot(2,2);...
        Trot(3,2)*Trot(3,3), Trot(3,3)*Trot(3,1), Trot(3,1)*Trot(3,2)];
T21=[Trot(2,1)*Trot(3,1), Trot(2,2)*Trot(3,2), Trot(2,3)*Trot(3,3);...
     Trot(3,1)*Trot(1,1), Trot(3,2)*Trot(1,2), Trot(3,3)*Trot(1,3);...
     Trot(1,1)*Trot(2,1), Trot(1,2)*Trot(2,2), Trot(1,3)*Trot(2,3)];
T22=[Trot(2,2)*Trot(3,3)+Trot(2,3)*Trot(3,2),...
     Trot(2,1)*Trot(3,3)+Trot(2,3)*Trot(3,1),...
     Trot(2,2)*Trot(3,1)+Trot(2,1)*Trot(3,2);...
     Trot(1,2)*Trot(3,3)+Trot(1,3)*Trot(3,2),...
     Trot(1,3)*Trot(3,1)+Trot(3,3)*Trot(1,1),...
     Trot(1,1)*Trot(3,2)+Trot(1,2)*Trot(3,1);...
     Trot(1,2)*Trot(2,3)+Trot(1,3)*Trot(2,2),...
     Trot(1,3)*Trot(2,1)+Trot(1,1)*Trot(2,3),...
     Trot(1,1)*Trot(2,2)+Trot(1,2)*Trot(2,1)];

% calculate the transpose of the overall rotation tensor
Ttran=[T11, T12; T21, T22];
% calculate the rotated elastic tensor
Cij=Ttran*Cij0*transpose(Ttran);

end