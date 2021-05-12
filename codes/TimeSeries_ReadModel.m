function rawmodel=TimeSeries_ReadModel(rawmodel,Nmodel)
%% Documentation
% Function read the input model
% Added by Xiaoran Chen on 10/23/2020
% Matlab R2016a

% If multiple model files exist, choose the one by specifying Nmodel.
% Naming convention of txt file is strict and it is has to be
% 'model-??.txt' where ?? stand for a number. 
filename=sprintf('model-%d.txt',Nmodel);
M=dlmread(filename,'',1,0);
fid = fopen(filename);
Parameterization = {fgetl(fid)};
fclose(fid);

if strcmp(Parameterization,'BackusNotation') % Backus Notation
    filename=sprintf('model-%d.txt',Nmodel);
    fid = fopen(filename);
    Model=textscan(fid,'%f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fid);
    rawmodel.z=flip(Model{1});
    rawmodel.alpha=flip(Model{2});
    rawmodel.beta=flip(Model{3});
    rawmodel.rho=flip(Model{4});
    rawmodel.xi_deg=flip(Model{5});
    rawmodel.tilt_deg=flip(Model{6});
    rawmodel.B_ani=flip(Model{7});
    rawmodel.E_ani=flip(Model{8});
    rawmodel.parameterization='BackusNotation';

elseif strcmp(Parameterization,'36-componentTensor') % Full tensor, 6 by 6 matrix
    filename=sprintf('model-%d.txt',Nmodel);
    M0=dlmread(filename,'',1,0);
    n=0;
    size_ans=size(M0);
    for i=(size_ans(1)-1)/7:-1:1
        rawmodel.z(n+1)=flip(M(7*n+1,1));
        rawmodel.alpha(n+1)=flip(M(7*n+1,2));
        rawmodel.beta(n+1)=flip(M(7*n+1,3));
        rawmodel.rho(n+1)=flip(M(7*n+1,4));
        rawmodel.ang1(n+1)=flip(M(7*n+1,5));
        rawmodel.ang2(n+1)=flip(M(7*n+1,6));
        rawmodel.ang3(n+1)=flip(M(7*n+1,7));
        rawmodel.C(:,:,i)=M((7*n+2):(7*n+7),1:6);
        n=n+1;
    end
    rawmodel.parameterization='36-componentTensor';

elseif strcmp(Parameterization,'81-componentTensor') % Full tensor, 9 by 9 matrix
    filename=sprintf('model-%d.txt',Nmodel);
    M0=dlmread(filename,'',1,0);
    n=0;
    size_ans=size(M0);
    for i=(size_ans(1)-1)/10:-1:1
        rawmodel.z(n+1)=flip(M0(10*n+1,1));
        rawmodel.alpha(n+1)=flip(M0(10*n+1,2));
        rawmodel.beta(n+1)=flip(M0(10*n+1,3));
        rawmodel.rho(n+1)=flip(M0(10*n+1,4));
        rawmodel.ang1(n+1)=flip(M0(10*n+1,5));
        rawmodel.ang2(n+1)=flip(M0(10*n+1,6));
        rawmodel.ang3(n+1)=flip(M0(10*n+1,7));

        C0=M0((10*n+2):(10*n+10),1:9);

        transA=reshape(C0',81,[]);
        a=zeros(3,3,3,3);
        m=1;

        for ll=1:3
            for kk=1:3
                for jj=1:3
                    for ii=1:3
                        a(ii,jj,kk,ll)=transA(m);
                        m = m + 1;
                    end
                end
            end
        end
        rawmodel.Cijkl(:,:,i)=C0;
        rawmodel.C(:,:,i)=Full2Voigt(a);
        n=n+1;
    end
    rawmodel.parameterization='81-componentTensor';

elseif strcmp(Parameterization,'ThomsenNotation') % Thomsen Parameters
    filename=sprintf('model-%d.txt',Nmodel);
    fid = fopen(filename);
    Model=textscan(fid,'%f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fid);
    rawmodel.z=flip(Model{1});
    rawmodel.alpha=flip(Model{2});
    rawmodel.beta=flip(Model{3});
    rawmodel.rho=flip(Model{4});
    rawmodel.epislon=flip(Model{5});
    rawmodel.gamma=flip(Model{6});
    rawmodel.delta=flip(Model{7});

    rawmodel.C=zeros(6,6,length(rawmodel.z));
    lambda=rawmodel.rho(1)*rawmodel.alpha(1)*rawmodel.alpha(1)-2*rawmodel.rho(1)*rawmodel.beta(1)*rawmodel.beta(1);
    mu=rawmodel.rho(1)*rawmodel.beta(1)*rawmodel.beta(1);
    lame=rawmodel.rho(1)*rawmodel.alpha(1)*rawmodel.alpha(1);
    rawmodel.C(:,:,1)=[lame lambda lambda 0 0 0;...
                lambda lame lambda 0 0 0;...
                lambda lambda lame 0 0 0;...
                 0 0 0 mu 0 0;...
                 0 0 0 0 mu 0;...
                 0 0 0 0 0 mu];
    for i=2:length(rawmodel.z)
        C33=rawmodel.rho(i)*rawmodel.alpha(i)*rawmodel.alpha(i);
        C44=rawmodel.rho(i)*rawmodel.beta(i)*rawmodel.beta(i);
        C11=(2*rawmodel.epislon(i)+1)*C33;
        C66=(2*rawmodel.gamma(i)+1)*C44;
        C13=sqrt((C33-C44)*(2*rawmodel.delta(i)*C33+C33-C44))-C44;

        rawmodel.C(:,:,i)=[C11 C11-2*C66 C13 0 0 0;...
                    C11-2*C66 C11 C13 0 0 0;...
                    C13 C13 C33 0 0 0;...
                     0 0 0 C44 0 0;...
                     0 0 0 0 C66 0;...
                     0 0 0 0 0 C66];
    end
    rawmodel.para_hs=[rawmodel.alpha(1) rawmodel.beta(1) rawmodel.rho(1)];
    rawmodel.parameterization='ThomsenNotation';
end