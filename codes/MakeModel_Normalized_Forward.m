function model=MakeModel_Normalized_Forward(rawmodel)
%% Documentation
% Function to normalize the parameters in the model 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
density_n=5515;
length_n=6371000;
frequency_n=1.075190645e-3;
velocity_n=frequency_n*length_n;
elastic_n=density_n*velocity_n*velocity_n;

model=struct;
if strcmp(rawmodel.parameterization,'BackusNotation')
    for i=1:length(rawmodel)

            mu=rawmodel(i).rho.*rawmodel(i).beta.*rawmodel(i).beta;%lame parameter II
            lambda=rawmodel(i).rho.*rawmodel(i).alpha.*rawmodel(i).alpha-2*mu;%lame parameter I

            AA=lambda+2*mu;
            DD=mu;
            z=(rawmodel(i).z-rawmodel(i).z(1)*ones(size(rawmodel(i).z)))./length_n;%depth of interface 
            rho=rawmodel(i).rho./density_n;%density

            AA=AA./elastic_n;
            BB=rawmodel(i).B_ani.*AA;
            CC=0*AA;% same order for the 4 lobe coefficients
            DD=DD./elastic_n;
            EE=rawmodel(i).E_ani.*DD;

            para=[AA'; BB'; CC'; DD'; EE'; rawmodel(i).xi_deg'; rho'; rawmodel(i).tilt_deg'];
            para_hs=[rawmodel(i).alpha(1) rawmodel(i).beta(1) rawmodel(i).rho(1)];

            model(i).AA=AA;
            model(i).BB=BB;
            model(i).CC=CC;
            model(i).DD=DD;
            model(i).EE=EE;
            model(i).xi_deg=rawmodel(i).xi_deg;
            model(i).rho=rho;
            model(i).tilt_deg=rawmodel(i).tilt_deg;
            model(i).para=para;
            model(i).para_hs=para_hs;
            model(i).z=z;

    end
%     correct for the intital pulse
        model(i).tsum=0;

        deltaz=zeros(1,length(rawmodel(i).z)-1);
        for j=1:length(deltaz)
            deltaz(j)=rawmodel(i).z(j)-rawmodel(i).z(j+1);
            if strcmp(flag,'P')
               model(i).tsum=model(i).tsum+deltaz(j)/rawmodel(i).alpha(j+1);
            else
               model(i).tsum=model(i).tsum+deltaz(j)/rawmodel(i).beta(j+1);
            end
        end
        
else
    for i=1:length(rawmodel)

            para_hs=[rawmodel(i).alpha(1)./velocity_n rawmodel(i).beta(1)./velocity_n rawmodel(i).rho(1)./density_n];

            model(i).rho=rawmodel(i).rho./density_n;
            model(i).para_hs=para_hs;
            model(i).z=(rawmodel(i).z-rawmodel(i).z(1)*ones(size(rawmodel(i).z)))./length_n;%depth of interface 
            model(i).C(:,:,:)=rawmodel.C(:,:,:)./elastic_n;

    end
end
model.parameterization = rawmodel.parameterization;
    