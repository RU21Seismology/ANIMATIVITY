function stup=TimeSeries_ReadSetup(stup)
%% Documentation
% Function to read the set up information
% Added by Xiaoran Chen on 10/23/2020
% Matlab R2016a

fid = fopen('setup.txt');
tline = fgetl(fid);
ischar(tline)
C=textscan(tline,'%s %s %s');

stup.targetdir=char(C{1});% target directory where all files are stored (full path is required)
stup.sitename=char(C{2});% site name
stup.inputformat=char(C{3}); % input format of the model ('txt' or 'mat')
fclose(fid);

