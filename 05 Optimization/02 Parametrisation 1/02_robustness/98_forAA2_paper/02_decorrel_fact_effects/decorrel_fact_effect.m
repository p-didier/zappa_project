clear; close all; clc;
addpath('..\..')
addpath(genpath('..\..\..\..\02 MATLAB\basefunctions'))
addpath(genpath('..\..\..\toolbox'))
addpath(genpath('..\..\..\..\02 MATLAB\09 Reverberation Time sims\10 RTana_v3\toolbox'))

% -- Purpose of script
% Properly showing the effect of the decorrelation factor.

% (c) Paul Didier - 23-Mar-2023 11:31
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

% Data location and referencing
prePath = '..\..\..\..\01 ANSYS\07 Optimization\02_parametrisation1\02_exports\backups';
runRef_all = {'kwwcwpj','snqcogt','iillxzc','tyickyz','drmmvgu'};
runRef_all = {'snqcogt','iillxzc','tyickyz','drmmvgu'};
runRef_all = {'kwwcwpj','rxztxfn'};

% RNG seeds
seeds = 1:100;

% Forced physical parameters (to not force, use NaN)
ds = 0.2;       % Sample thickness

rng('default')   % SET INITIAL RANDOM GENERATOR SEED

%% READ DATA

combinedDataForPost = cell(length(runRef_all));
for idxmain = 1:length(runRef_all)
    runRef = runRef_all{idxmain};

    % Find appropriate folder
    if length(runRef) == 7
        folders = dir([prePath '\series_' runRef]);
        folders(1:2) = []; folders(~cell2mat({folders.isdir})) = [];
        designDataPath = [prePath '\series_' runRef '\' folders(end).name];
    else
        designDataPath = [prePath '\' runRef];
    end

    % Read data
    designData = readdesigndata(designDataPath,'Workspace',1,'OptiOutcome',1);

    % Edit/add forced parameters
    if exist('ds','var')
        if ~isnan(ds)
            designData.params.ds = ds;
        end
    end

    % Run
    status = run_with_various_seeds(designData,seeds);

    save(['seedstest_' runRef], 'status')
    
end
disp('ALL DONE.');

