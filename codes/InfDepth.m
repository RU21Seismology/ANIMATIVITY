function D0=InfDepth(Slowness,z, w0)
%% Documentation
% Function to calculate the initial depth dependent exponential matrix 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
    global Num
    
    D0=zeros(6,6,Num);
    
    %calculate the phase and amplitude spectrum for each frequency
    for cnt=1:Num
       
        %calculate the depth dependent matrix and depth difference       
        m0=-1i*w0(cnt)*z;
        n0=[exp(m0*Slowness(1)) exp(m0*Slowness(2)) exp(m0*Slowness(3)) exp(m0*Slowness(4)) exp(m0*Slowness(5)) exp(m0*Slowness(6))];
        D0(:,:,cnt)=diag(n0);
                
    end
end
