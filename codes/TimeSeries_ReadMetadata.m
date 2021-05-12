function info=TimeSeries_ReadMetadata(info,path)
%% Documentation
% Function to read the metadata
% Added by Xiaoran Chen on 10/23/2020
% Matlab R2016a

fid = fopen([path '/info.txt']);
tline = fgetl(fid);
i=1;
while ischar(tline)
    C=textscan(tline,'%s %f %f %f %f %f %d %f %s %f %f');
    info(i).SiteName=C{1};% site name
    info(i).dep_km=C{2};% depth of the source in kilometers
    info(i).dist_deg=C{3};% distance of the source in degrees
    info(i).baz_deg=C{4};% back azimuth of the source in degrees
    info(i).delta_sec=C{5};% sampling rate in seconds
    info(i).dur_sec=C{6};% duration of the record in seconds
    info(i).source=C{7};% source function expressed in its mathematical format
    info(i).wvlen_sec=C{8};% wave length in seconds
    info(i).flag=C{9};% type of the incident wave (1:P; 2:SV; 3:SH)
    info(i).lenplt_sec=C{10};% duration of the time series to plot
    info(i).pole=C{11};%the pole where the first pulse should be around
    tline = fgetl(fid);
    i=i+1;
end
fclose(fid);

