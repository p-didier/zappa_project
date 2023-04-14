clear; close all; clc;

% -- Purpose of script
% Manually edit the output file "RC_rxztxfn_SampleDims.mat" after having
% received the new data from Cédric.

% (c) Paul Didier - 28-Feb-2023 12:52
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

matFilename = '01_mat\99_backups\backup_20230228\RC_rxztxfn_SampleDims__old.mat';
pathToRefData = 'C:\Users\pdidier\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\07 Data files\00_detSEA';
destinationPath = '01_mat/RC_rxztxfn_SampleDims.mat';

load(matFilename)

%% EDIT STRUCTURE

newStruct = status;  % init
for ii = 1:size(newStruct.alphab, 2)
    % Get dims
    currDims = newStruct.RCparam(ii, :);
    % Reference data file name
    refDataFilename = ['dims' num2str(round(currDims(1) * 100)) 'x' ...
        num2str(round(currDims(2) * 100)) 'x20_Xi11k.m'];
    % Get ref data
    run([pathToRefData '/' refDataFilename])
    newStruct.alpharef(:, ii) = alpharef(1:size(newStruct.alpharef, 1));
    stop = 1;
end

%% SAVE TO MAT

status = newStruct;
save(destinationPath, 'check', 'runRef', 'status')

%% PLOT

close all
fig = figure; fig.Units = "Normalized"; fig.Position = [0.1 0.1 0.7 0.4];
subplot(121)
hold on; grid on
plot(newStruct.alpharef)
ylim([0, 1.4])
title 'new data'
subplot(122)
hold on; grid on
plot(status.alpharef)
ylim([0, 1.4])
title 'old data'

%% FUNCTIONS

