clear; close all; clc;

% -- Purpose of script
% Post-process simulation data for samples-dimensions robustness
% check (ZAPPA AA2 paper).

% (c) Paul Didier - 12-Mar-2023 16:58
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

runRefs = {'kwwcwpj','snqcogt','iillxzc','tyickyz','drmmvgu','rxztxfn'};

pathToData = '01_mat';

debugBool = 1;  % if 1, debug mode (20.03.2023)

%% LOAD DATA

allData = cell(size(runRefs));
for ii = 1:length(runRefs)
    currFilename = [pathToData '/RC_' runRefs{ii} '_SampleDims.mat'];
    load(currFilename)
    allData{ii} = status;
end

%% PLOT

if debugBool
    fig = figure; fig.Units = "Normalized"; fig.Position = [0.1 0.1 0.7 0.8];
    for ii = 1:length(allData)
        subplot(2, 3, ii)
        hold on; grid on
        plot(1:4, allData{ii}.alphab)
        plot(1:4, allData{ii}.alpharef, 'k')
%         delta = abs(allData{ii}.alphab(:, 2:end) - allData{ii}.alpharef(:, 2:end));
%         boxplot(delta.','Labels',{'100', '125', '160', '200'})
%         scatter(1:4, abs(allData{ii}.alphab(:, 1) - allData{ii}.alpharef(:, 1)), 10, 'k')
        xlim([0.5, 4.5])
        ylim([0, Inf])
        xlabel '$f_\mathrm{c}$ [Hz]'
        ylabel '$\Delta_\mathrm{a}\alpha$'
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        ylabel '$\Delta_\mathrm{a}\alpha$'
        title(['Opt. design ' num2str(ii)])
        if ii ~= 1 && ii ~= 4
            ylabel('')
        end
        ylim([0.85, 1.3])
    end
else
    fig = figure; fig.Units = "Normalized"; fig.Position = [0.1 0.1 0.7 0.8];
    for ii = 1:length(allData)
        subplot(2, 3, ii)
        hold on; grid on
        delta = abs(allData{ii}.alphab(:, 2:end) - allData{ii}.alpharef(:, 2:end));
        boxplot(delta.','Labels',{'100', '125', '160', '200'})
        scatter(1:4, abs(allData{ii}.alphab(:, 1) - allData{ii}.alpharef(:, 1)), 10, 'k')
        xlim([0.5, 4.5])
        ylim([0, Inf])
        xlabel '$f_\mathrm{c}$ [Hz]'
        ylabel '$\Delta_\mathrm{a}\alpha$'
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        ylabel '$\Delta_\mathrm{a}\alpha$'
        title(['Opt. design ' num2str(ii)])
    %     if ii <= 3
    %         xlabel('')
    %         xticklabels({})
    %     end
        if ii ~= 1 && ii ~= 4
            ylabel('')
    %         yticklabels({})
        end
        ylim([0, 0.35])
    %     set(gca, 'PlotBoxAspectRatio', [1,1,1])
    end
end

if 1
    if debugBool
        exportfigure(gcf, 'sampleDimsRob_debug', '00_fig/03_sampleDimsAll')
    else
        exportfigure(gcf, 'sampleDimsRobustness', '00_fig/03_sampleDimsAll')
    end
end