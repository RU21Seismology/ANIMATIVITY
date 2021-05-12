function V0=Spec0(Amp,initime,w)
%% Documentation
% Function to initialize the spetrum of the amplitude/constant part of the matrix 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here

    global Num length_n frequency_n
    V0=zeros(6,6,Num);    

    A1=ones(1,Num)./length_n*Amp(1);
    P1=-initime(1)*frequency_n.*w;    
    
    A2=ones(1,Num)./length_n*Amp(2);
    P2=-initime(2)*frequency_n.*w;    

    A3=ones(1,Num)./length_n*Amp(3);
    P3=-initime(3)*frequency_n.*w;    

    %calculate the initial spectrum using amplitude and phase spectrum
    V0(4,4,:)=A1.*(exp(P1*1i));
    V0(5,5,:)=A2.*(exp(P2*1i));
    V0(6,6,:)=A3.*(exp(P3*1i));
    
end
   
   