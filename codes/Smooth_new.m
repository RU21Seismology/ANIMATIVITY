function B=Smooth_new(B_old,df,sps)
%% Documentation
% Function to multiply all the results with a half the period of cos.^2
% function to make it smoother 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
global Num frequency_n

B=zeros(6,6,Num);
%use 50 sps data so that the Nyquist is 25Hz
fmax=sps/2/frequency_n;

dw=2*pi*df;

for cnt=1:Num
  w=dw*(cnt-1);
  if w>(fmax*2*pi)
    B(:,:,cnt)=0;
  else
    B(:,:,cnt)=B_old(:,:,cnt)*(1+cos(w/(2*fmax)))*(1+cos(w/(2*fmax)))/4;
  if cnt == 1
    B(:,:,cnt)=0;
  end
  end
end


