function wvfm=PostProcess_TimeSeries(info,model,stup,rawmodel,wvfm,Nmodel)
%% Documentation
% Function to generate the synthetics and save the ourput
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here

for j=1:length(info)

    [R,T,Z,raypar]=SynthWF(info(j),model);
    wvfm(j).SiteName=info(j).SiteName;
    wvfm(j).dep_km=info(j).dep_km;
    wvfm(j).dist_deg=info(j).dist_deg;
    wvfm(j).baz_deg=info(j).baz_deg;
    wvfm(j).delta_sec=info(j).delta_sec;
    wvfm(j).dur_sec=info(j).dur_sec;
    wvfm(j).pole=info(j).pole;
        
    if info(j).source == 1
        wvfm(j).sourcefunc='sin^2';
    elseif info(j).source == 2
        wvfm(j).sourcefunc='sincos';
    elseif info(j).source == 3
        wvfm(j).sourcefunc='rand1';
    elseif info(j).source == 4
        wvfm(j).sourcefunc='rand2';
    end
  
    wvfm(j).phase = info(j).flag;

    dur=2^nextpow2(wvfm(j).dur_sec/wvfm(j).delta_sec);
    temp=info(j).delta_sec*(0:dur-1);
    wvfm(j).time=temp';
    K= ceil(info(j).lenplt_sec/wvfm(j).delta_sec)+1;
    wvfm(j).wvlen_sec=info(j).wvlen_sec;

    wvfm(j).model=rawmodel;
    wvfm(j).R=R';
    wvfm(j).T=T';
    wvfm(j).Z=Z';
   
    waveform = wvfm(j);
    Savepdf(wvfm(j).time,R,T,Z,wvfm(j).baz_deg,K,char(waveform.SiteName),Nmodel,stup.targetdir);

    wvfm(j).Rname=sprintf('%s-%d-%.2f-%.1f-R.txt',char(waveform.SiteName), Nmodel, raypar, wvfm.baz_deg);
    wvfm(j).Tname=sprintf('%s-%d-%.2f-%.1f-T.txt',char(waveform.SiteName), Nmodel, raypar, wvfm.baz_deg);
    wvfm(j).Zname=sprintf('%s-%d-%.2f-%.1f-Z.txt',char(waveform.SiteName), Nmodel, raypar, wvfm.baz_deg);

end