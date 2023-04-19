clear; close all; clc;

% -- Purpose of script
% Post-process data from decorrelation factor effect test.

% (c) Paul Didier - 24-Mar-2023 08:27
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

testType = 'changing_seed';  % repeating the same simulation with a 
    % different rng(seed) call at every run.
testType = 'changing_state'; % repeating the same simulation with 
    % the same rng(seed) (outside of the for-loop), but rand() is
    % called at every run --> the rng state changes.

runRefs = {'kwwcwpj','snqcogt','iillxzc','tyickyz','drmmvgu','rxztxfn'};
% runRefs = {'snqcogt','iillxzc','tyickyz','drmmvgu'};
matFilenamePrefix = ['./out/' testType '_'];
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
        legend([l1(1), l2, l3], {'$\alpha_\mathrm{sim}$ for individual runs',...
            'Mean and standard dev. of $\alpha_\mathrm{sim}$ across runs',...
            '$\alpha_\mathrm{diff}$ (reference)'},...
            'Location','southeast')
    end
    OTOBticks
    ylabel '$\alpha$'
    ylim([0.9, 1.3])
    title(['RR\#' num2str(ii)])
end
if strcmp(testType, 'changing_seed')
    testTypeStr = 'Changing rng(seed)';
elseif strcmp(testType, 'changing_state')
    testTypeStr = 'Changing state, with rng(0)';
end
suptitle([testTypeStr ' -- ' num2str(nSeedsToConsider) ' runs'])

if 0
    exportfigure(fig, testType, 'fig')
end
