clear; close all; clc;

% -- Purpose of script
% Adapt .mat file for design "rxztxfn" --> discard [2.83, 4.00] m^2 sample
% and add [2.90, 4.10] m^2 sample.

% (c) Paul Didier - 02-Mar-2023 11:02
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

% Load new results
load results_newsample
% Load old .mat file
pathToMatFile = '..\..\01_mat\99_backups\backup_20230302\RC_rxztxfn_SampleDims.mat';
load(pathToMatFile)
% Export path
exportPath = '..\..\01_mat\RC_rxztxfn_SampleDims.mat';

idxFaultyDims = 12;
newDims = [2.90039576600427,4.10304146667245];

%% PROCESS

status_new = status;  % copy
status_new.alphab(:, idxFaultyDims) = results.alphasimb;
status_new.alpharef(:, idxFaultyDims) = results.alpharef;
status_new.RCparam(idxFaultyDims, :) = newDims;

%% SAVE

status = status_new;
save(exportPath, 'check', 'runRef', 'status')
