function [R,T,Z]=GetDs(B)
%% Documentation
% Function to extract the displacement from the displacement-stress matrix 
% and transform back into time domain
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here    
    global length_n
    %downgoing P wave
    ux_dp=squeeze(B(1,1,:));
    uy_dp=squeeze(B(2,1,:));
    uz_dp=squeeze(B(3,1,:));
    
    %downgoing S1 wave
    ux_ds1=squeeze(B(1,2,:));
    uy_ds1=squeeze(B(2,2,:));
    uz_ds1=squeeze(B(3,2,:));
    
    %downgoing S2 wave
    ux_ds2=squeeze(B(1,3,:));
    uy_ds2=squeeze(B(2,3,:));
    uz_ds2=squeeze(B(3,3,:));
    
    %upgoing P wave
    ux_up=squeeze(B(1,4,:));
    uy_up=squeeze(B(2,4,:));
    uz_up=squeeze(B(3,4,:));
    
    %upgoing S1 wave
    ux_us1=squeeze(B(1,5,:));   
    uy_us1=squeeze(B(2,5,:)); 
    uz_us1=squeeze(B(3,5,:));
    
    %upgoing S2 wave
    ux_us2=squeeze(B(1,6,:));
    uy_us2=squeeze(B(2,6,:));
    uz_us2=squeeze(B(3,6,:));
   
    %downgoing P wave
    ux_dp=ifft(ux_dp,'symmetric').*length_n;
    uy_dp=ifft(uy_dp,'symmetric').*length_n;
    uz_dp=ifft(uz_dp,'symmetric').*length_n;
    
    %downgoing S1 wave
    ux_ds1=ifft(ux_ds1,'symmetric').*length_n;
    uy_ds1=ifft(uy_ds1,'symmetric').*length_n;
    uz_ds1=ifft(uz_ds1,'symmetric').*length_n;
    
    %downgoing S2 wave
    ux_ds2=ifft(ux_ds2,'symmetric').*length_n;
    uy_ds2=ifft(uy_ds2,'symmetric').*length_n;
    uz_ds2=ifft(uz_ds2,'symmetric').*length_n;
    
    %upgoing P wave
    ux_up=ifft(ux_up,'symmetric').*length_n;
    uy_up=ifft(uy_up,'symmetric').*length_n;
    uz_up=ifft(uz_up,'symmetric').*length_n;
    
    %upgoing S1 wave
    ux_us1=ifft(ux_us1,'symmetric').*length_n;
    uy_us1=ifft(uy_us1,'symmetric').*length_n;
    uz_us1=ifft(uz_us1,'symmetric').*length_n;
    
    %upgoing S2 wave
    ux_us2=ifft(ux_us2,'symmetric').*length_n;
    uy_us2=ifft(uy_us2,'symmetric').*length_n;
    uz_us2=ifft(uz_us2,'symmetric').*length_n;

    R=ux_dp+ux_ds1+ux_ds2+ux_up+ux_us1+ux_us2;
    antiT=uy_dp+uy_ds1+uy_ds2+uy_up+uy_us1+uy_us2;
    antiZ=uz_dp+uz_ds1+uz_ds2+uz_up+uz_us1+uz_us2;

    T=-antiT;
    Z=-antiZ;
    
    R=R.';
    T=T.';
    Z=Z.';

end
 

