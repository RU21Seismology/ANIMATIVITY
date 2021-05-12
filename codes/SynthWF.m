function [Rcomp,Tcomp,Zcomp,raypar]=SynthWF(info,model)
%% Documentation
% Function to calculate the 3-component synthetic waveforms
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Provide Normalization coefficients
% This set of normal mode normalization coefficients are inherited from
% anirec.f by Jeffrey Park and Vadim Levin to reduce the condition numbers.
% Tests on a couple of parameters show that even without the
% normalizations, the code will run stably but it won't hurt to keep it.
% Units here are SI units. 
global frequency_n length_n elastic_n density_n velocity_n
density_n=5515;
length_n=6371000;
frequency_n=1.075190645e-3;
velocity_n=frequency_n*length_n;
elastic_n=density_n*velocity_n*velocity_n;

dep=info.dep_km;
gcarc=info.dist_deg;
baz=info.baz_deg;
dt=info.delta_sec;
flag=info.flag;
source=info.source;
per=info.wvlen_sec;
T=info.dur_sec;

%% Set up the intervals and durations in time and frequency domains
global p;
if strcmp(flag,'P')
   tt=taupTime('iasp91',dep,'P,PKP,PKIKP,PKiKP','deg',gcarc);
elseif strcmp(flag, 'SV')
   tt=taupTime('iasp91',dep,'S,SKS','deg',gcarc);
elseif strcmp(flag, 'SH')
   tt=taupTime('iasp91',dep,'S,SKS','deg',gcarc);
end

found = 0;
    
if found == 0
    for i=1:length(tt)
        if strcmp(tt(i).phaseName,'P')
            ttp = tt(i);
            found = 1;
            break
        end
    end
end

if found == 0
    for i=1:length(tt)
        if strcmp(tt(i).phaseName,'PKP')
            ttp = tt(i);
            found = 1;
            break
        end
    end
end

if found == 0
    for i=1:length(tt)
        if strcmp(tt(i).phaseName,'PKIKP')
            ttp = tt(i);
            found = 1;
            break
        end
    end
end
    
if found == 0
    for i=1:length(tt)
        if strcmp(tt(i).phaseName,'PKiKP')
            ttp = tt(i);
            found = 1;
            break
        end
    end
end

if found == 0
    disp('Phase not found');
end

r2d=pi/180;
raypar=ttp.rayParam*r2d;
p=raypar/111000*velocity_n;
sps=1/dt;

% normally a global variable may be confusing and we do not recommand for
% most of the time, but here I find it pretty useful
%number of points in calculation
global Num;
num=T*sps;
%for the Fourier transform we need num=pow(2,n)
Num=2^nextpow2(num);

% Preallocate the memory for different array to bring up the speed
% Displacement stress matrix for six phases
B_before=zeros(6,6,Num);
% Displacement stress matrix convolved with a source
B_before_conv=zeros(6,6,Num);

%Calculate the frequency spacing and normalize it
f_begin=0;
%interval in frequency domain is the reciptical of duration of the records
df=1/((Num-1)*dt)/frequency_n;
w=2*pi*(f_begin+df*(0:Num-1));

if strcmp(model.parameterization,'BackusNotation')
    z=model.z;
    baz_deg=baz*ones(1,length(z));% from N
    model.baz_deg=baz_deg;
    para=[model.para; model.baz_deg];
    para_hs=model.para_hs;
    rho=model.rho;
    Cij=BuildGeneral(para);
else
    z=model.z;
    para_hs=model.para_hs;
    rho=model.rho;
    Cij=model.C;
end
% number of layers including both the top layer and the half space
N=length(z);
% each element is the vertical slowness of six phases
Slowness=zeros(1,6,N);
% each element is the polarization of six phases/norm of the phase front
Polarization=zeros(3,6,N);

% calculate the phase velocities and polarization vectors for both the half
% space and the top layer
[Slowness_hs, Polarization_hs]=EigHs(para_hs);
[Slowness_top, Polarization_top]=QEP_general(Cij(:,:,end),rho(end));

% you can move the three following lines to the beginning of the function.

if strcmp(flag,'P')
    Amp=[1 0 0];
    initime=[info.pole - model.tsum 0 0];
elseif strcmp(flag,'SV')
    Amp=[0 1 0];
    initime=[0 info.pole - model.tsum 0];
elseif strcmp(flag, 'SH')
    Amp=[0 0 1];
    initime=[0 0 info.pole - model.tsum];    
end

% calculate the initial spectrum of the input wave (product of the amplitude spectrum and the phase spectrum)
V0=Spec0(Amp,initime,w); 
D0=InfDepth(Slowness_hs,z(1),w);

% since the hs (half space) involves a slightly differnt calculation
% routine, hereafter we separate it from the rest of the computation
Slowness(:,:,1)=Slowness_hs;
Polarization(:,:,1)=Polarization_hs;

for i=2:N
    [Slowness(:,:,i), Polarization(:,:,i)]=QEP_general(Cij(:,:,i),rho(i));
end

%here we start to loop over frequency
% for cnt=1:2000%check if we need such a large number
for cnt=1:Num/2
    % compute the eigenvalue matrix for displacement and stress
    E_top=EigStr_general(Slowness_top, Polarization_top, Cij(:,:,end), w(cnt)); % top, free surface

    % compute the initial spetrum and make it diagonal
    W_hs=diag(D0(:,:,cnt)*V0(:,:,cnt));
    Wu_hs=W_hs(4:6,:); % only use the upgoing wave (Wu)
    
    if w(cnt)==0%try without zeros
       Wu=Wu_hs;
       Wd=[0;0;0];
       deltaD=diag([1 1 1 1 1 1]);
    else
       % for non-zero phase, equations can be found in Chen, 1993; Park and
       % Levin, 1998; Park and Levin, 2016...
       % Generalized reflection coefficients from upgoing waves
       Rud0=Ds_Inf_FS(E_top, w(cnt), z(end-1), z(end), Slowness_top);
       RUD0_hat_fs=Rud0;
    
       for i=N-1:-1:1
          % eigenvalue matrix for displacement-stress in upper and lower
          % layers
          E2=EigStr_general(Slowness(:,:,i), Polarization(:,:,i), Cij(:,:,i), w(cnt));
          E1=EigStr_general(Slowness(:,:,i+1), Polarization(:,:,i+1), Cij(:,:,i+1), w(cnt));

          if i==1
             % compute the generalized R/T coefficients generated by
             % upgoing waves
             [TU1, RUD1]=InfRecursive(E1, E2, w(cnt), z(i+1), z(i), z(i), Slowness(:,:,i+1), Slowness(:,:,i), Rud0);

          else
             [TU1, RUD1]=InfRecursive(E1, E2, w(cnt), z(i+1), z(i-1), z(i), Slowness(:,:,i+1), Slowness(:,:,i), Rud0);          
          end
          TU_hat(:,:,i)=TU1;
          Rud0=RUD1;% 0 stands for free surface
          
       end

      Wu=Wu_hs;
      for i=1:N-1
         % calculate the upgoing waves approaching the free surface
         Wu=TU_hat(:,:,i)*Wu;
      end
      % downgoing waves generated at the free surface from the upgoing waves
      Wd=RUD0_hat_fs*Wu;
      
      z_expect=z(end);
      z0=z(end);
      z1=z(end-1);

      delta1_d_z0=[exp(-1i*w(cnt)*Slowness_top(1)*(z_expect-z0)) exp(-1i*w(cnt)*Slowness_top(2)*(z_expect-z0)) exp(-1i*w(cnt)*Slowness_top(3)*(z_expect-z0))];
      delta1_u_z0=[exp(-1i*w(cnt)*Slowness_top(4)*(z_expect-z1)) exp(-1i*w(cnt)*Slowness_top(5)*(z_expect-z1)) exp(-1i*w(cnt)*Slowness_top(6)*(z_expect-z1))];

      %calculate the phase changes according to the depth
      deltaD=diag([delta1_d_z0 delta1_u_z0]);
    end
      %compose the displacement-stress matrix after wave propagation
      B_before(:,:,cnt)=E_top*deltaD*diag([Wd.' Wu.']);

end

% generate a sin^2 function to simulate the source
ImpulseFun=zeros(1,Num);
if source==1
    for ii=1:ceil(per/dt)
        tt=dt*(ii-1);
        ImpulseFun(ii)=sin(pi/per*tt)^2;
    end
% generate a sincos function to simulate the source
elseif source==2
    for ii=1:ceil(per/dt)
        tt=dt*(ii-1);
        ImpulseFun(ii)=sin(pi/per*tt)*cos(pi/per*tt);
    end
% generate a 10s random-walk function to simulate the source
elseif source==3
    per = 10;
    nwa = ceil(per/dt);
    xx = 0;
    frrr = 1;
    for ii=1:nwa
        xx = xx + rand*2*dt;
        phase = 2 *pi * frrr*xx;
        ImpulseFun(ii)=(sin(phase)+sin(2.66*phase))*exp(-0.02*xx^2);
    end
% generate a 10s random walk function to simulate the source
elseif source==4
    per = 10;
    nwa = ceil(per/dt);
    xx = 0;

    for ii=1:nwa
        xx = xx + rand-0.5;
        phase = pi/2*(ii-1)/nwa;
        ImpulseFun(ii)=cos(phase)*xx;
    end
end

% compute the spectrum of the source
IMPULSE=fft(ImpulseFun);

for i=1:6
    for j=1:6
        for cnt=1:Num
            % conv in time domain = product of spectrum
            B_before_conv(i,j,cnt)=B_before(i,j,cnt)*IMPULSE(cnt);
        end
    end
end

%% Smooth
% taper the spectrum to get rid of the ringing effects
B=Smooth_new(B_before_conv,df,sps);
% inverse-fourier transform to get the RTZ components
[Rcomp,Tcomp,Zcomp]=GetDs(B);



