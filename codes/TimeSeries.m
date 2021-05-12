% This is the main function to drive all the subroutines in ANIMATIVITY.

% ANIMATIVITY stands for 'ANIsotropic MATlab ReflectIVITY'. It computes the
% synthetic time series after upgoing seismic body waves propagate through
% the horizontally stratified anisotropic media. It is written based on the
% reflectivity algorithm proposed by Fuchs and Müller (1971). 

% ANIMATIVITY can accept four different descriptions of anisotropy:
% 1) Backus Notation (Backus, 1965)
% 2) Full description of the tensor (36 components)
% 3) Full description of the tensor (81 components)
% 4) Thomsen parameters (Thomsen, 1986)

% An earlier version of this code in Fortran named as 'anirec' (Levin and
% Park, 1997) focuses on Backus Notation and can be accessed from IRIS
% website http://seiscode.iris.washington.edu/projects/rfsyn.

% This software was developed by Xiaoran Chen with great help from Dr.
% Jeffrey Park and Dr.Vadim Levin.

% This software has been tested on Matlab 2016a.
% Updated on 10/23/2020.

clear;
clc;
close;

%% Read the information of environment setup including the I/O directory,
% site name, format of model input and incident wave type.
stup=struct;
stup=TimeSeries_ReadSetup(stup);

if ~exist([stup.targetdir '/Output'], 'dir')
    mkdir([stup.targetdir '/Output']);
end    

%% Read the metadata
info=struct;
info=TimeSeries_ReadMetadata(info,[stup.targetdir '/Input']);

Nmodel=0; % This number marks which model will be used in computation

if strcmp(stup.inputformat,'txt')   % read model from .txt file
    rawmodel=struct;
    rawmodel=TimeSeries_ReadModel(rawmodel,Nmodel);
    save([stup.targetdir '/Output/rawmodel.mat'],'rawmodel');
else % read model from .mat file
    rawmodel=load([stup.targetdir '/Input/rawmodel.mat']);
end

% normalize the model for the code input
model=MakeModel_Normalized_Forward(rawmodel);

%% Conduct forward modelling
wvfm=struct;
path=[stup.targetdir '/Output'];

% generate synthetics
wvfm=PostProcess_TimeSeries(info,model,stup,rawmodel,wvfm,Nmodel);
save([stup.targetdir '/Output/waveform.mat'],'wvfm');

%% Save corresponding files
% save the selected components of the synthetics with time axis.
dlmwrite([stup.targetdir '/Output/' sprintf('%s', wvfm(1).Rname)],wvfm(1).R);
