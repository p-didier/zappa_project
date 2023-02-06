clear; close all; clc;
addpath(genpath('..\..\..\..\..\02 MATLAB\basefunctions'))

% -- Purpose of script
% Figure edit for Applied Acoustics submission #2.

% (c) Paul Didier - 24-Nov-2022 11:53
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

figName = 'out_V250rxy85rxz140';
% figName = 'out_V250rxy114rxz190';

%% PROCESS

fig = openfig(figName);
ax = fig.Children(2);
fig.Position(4) = fig.Position(4) + 0.01;

% Edit x-axis label
ax.XLabel.String = '$\sigma$ [kNsm\textsuperscript{-4}]';

% Add lines
hold on
xLines = linspace(ax.XLim(1), 48700 + 2.5409e+03, size(ax.Children(1).CData, 2) + 1);
lw = 0.5;
xline(xLines, LineWidth=lw, Alpha=1, Color='k')
yline(normToExactOTOBs([100, 125, 160, 200]) * 2^(1/6), LineWidth=lw, Alpha=1, Color='k')

%% EXPORT

if 1
    exportfigure(gcf, [figName '_forAA2'])
end