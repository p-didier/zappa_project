clear; close all; clc;

% -- Purpose of script
% Post-process data from decorrelation factor effect test.

% (c) Paul Didier - 24-Mar-2023 08:27
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

runRefs = {'kwwcwpj','snqcogt','iillxzc','tyickyz','drmmvgu','rxztxfn'};
% runRefs = {'snqcogt','iillxzc','tyickyz','drmmvgu'};
matFilenamePrefix = 'seedstest_';
nSeedsToConsider = 50;

%% PROCESS

f = [100, 125, 160, 200];

close all
fig = figure; fig.Units = "Normalized"; fig.Position = [0.1 0.1 0.9 0.8];
for ii = 1:length(runRefs)
    matFilename = [matFilenamePrefix runRefs{ii} '.mat'];
    load(matFilename)
    
%     nSeeds = size(alpharef_all, 2);

    % Plot
    subplot(2,3,ii)
    hold on; grid on
    l1 = plot(f, status.alphab(:, 1:nSeedsToConsider), 'o-', 'Color', 0.8*ones(3,1));
    l2 = errorbar(f, mean(status.alphab(:, 1:nSeedsToConsider), 2), ...
        std(status.alphab(:, 1:nSeedsToConsider), [], 2), 'bo-', "linewidth", 1.5);
    l3 = plot(f, status.alpharef(:, 1), 'xk-', "linewidth", 1.5);
    if ii == length(runRefs)
        legend([l1(1), l2, l3], {'$\alpha_\mathrm{sim}$ for individual \texttt{rng} seeds',...
            'Mean and standard dev. of $\alpha_\mathrm{sim}$ across \texttt{rng} seeds',...
            '$\alpha_\mathrm{diff}$ (reference)'},...
            'Location','southeast')
    end
    OTOBticks
    ylabel '$\alpha$'
    ylim([0.9, 1.3])
    title(['RR\#' num2str(ii) ' - ref. "' runRefs{ii} '"'])
end
% suptitle(['Number of \texttt{rng} seeds considered: ' ...
%     num2str(nSeedsToConsider)])

if 1
    exportfigure(fig, 'rngTest', '00_fig/04_debugging')
end
