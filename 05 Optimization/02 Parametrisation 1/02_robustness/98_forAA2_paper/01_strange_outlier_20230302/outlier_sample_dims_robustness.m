clear; close all; clc;
addpath(genpath('..\..'))
addpath('..\..\..')
addpath(genpath('..\..\..\..\toolbox'))
addpath(genpath('..\..\..\..\..\02 MATLAB\09 Reverberation Time sims\10 RTana_v3\toolbox'))

% -- Purpose of script
% Single run of the [2.83, 4.00] m^2 sample (outlier -- see figure in
% SOUNDS PhD journal 2023, week 09, THU).

% (c) Paul Didier - 02-Mar-2023 10:25
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

runRef = 'rxztxfn';
check = 'SampleDims';
prePath = '..\..\..\..\..\01 ANSYS\07 Optimization\02_parametrisation1\02_exports\backups';


%% PROCESS

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

status = perf_robust_checks(designData,'Check',check);

%% PLOT

% close all
% fig = figure; fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.3646 0.4861];
% hold on; grid on

% ax = gca;

%% FUNCTIONS

