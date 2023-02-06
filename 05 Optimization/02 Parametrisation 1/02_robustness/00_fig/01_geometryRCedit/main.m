clear; close all; clc,
addpath(genpath('..\..\..\..\..\02 MATLAB\basefunctions'))

% -- Purpose of script
% Edit for the geometrical inaccuracies robustness check -- paper #1
% Internoise 2021.

% (c) Paul Didier - 13-May-2021 17:00

%% INIT

% runref = 'hdkrkyy';
% runref = 'wjiesbx';
% runref = 'fecjeep';
% runref = 'wxrayqc';
% runref = 'snqcogt';
% runref = 'iillxzc';
runref = 'tyickyz';
runref = 'kwwcwpj';
runref = 'drmmvgu';
runref = 'rxztxfn';

%% FETCH REFERENCE DATA

fnameFig = ['..\RC_' runref '_Geometry.fig'];
optiFig = ['..\..\..\00_fig\finalState_' runref '.fig'];

fig = openfig(optiFig,'invisible');
ax = fig.Children;
ax = ax(3);
data = ax.Children;
data = data(2);
alphaopti = data.YData;
fcopti = exactToNormOTOBs(data.XData); clear ax data
close(fig);

%% OPEN FIGURE

fig = openfig(fnameFig);
fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.45 0.2];
ax = fig.Children;
hold(ax,'on')
lines = ax.Children;

plot(ax,fcopti,alphaopti,'ro','linewidth',1)
ph = plotstairs(fcopti,alphaopti); ph.Color = 'red'; ph.LineWidth = 1;
ax.YLim = [0.75 1.35];
uistack(lines(1),'top')
ax.XLabel.String = '$f_\mathrm{c}$ [Hz]';
ax.Children(4).LineStyle = '-.';
ax.Children(4).LineWidth = ax.Children(2).LineWidth;
ax.Children(5).LineWidth = ax.Children(2).LineWidth;
% Delete outlier markers and median lines
idx = contains({lines(1).Children.Tag},{'Outliers','Median'});
delete(lines(1).Children(idx))   
% Turn whisker lines blue
idx = contains({lines(1).Children.Tag},{'Lower','Upper'});
idxt = 1:length(lines(1).Children);
idxt(~idx) = []; 
for ii = 1:length(idxt)
    lines(1).Children(idxt(ii)).Color = 'blue';
    lines(1).Children(idxt(ii)).LineStyle = '-';
end
ax.TickLabelInterpreter = 'latex';
%% EXPORT

if 0
    exportfigure(gcf, ['out_' runref], '00_fig')
end