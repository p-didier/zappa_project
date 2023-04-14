clear; close all; clc
addpath('..')
addpath(genpath('..\..\..\02 MATLAB\basefunctions'))
addpath(genpath('..\..\toolbox'))
addpath(genpath('..\..\..\02 MATLAB\09 Reverberation Time sims\10 RTana_v3\toolbox'))
% rmpath(genpath('..\..\_beforeFERR'))
% -- Purpose of script
% Assesses the robustness of a given RR design to various changes in
% experimental conditions. 

% (c) Paul Didier - 11-May-2021 10:16

%% INIT

% Data location and referencing
prePath = '..\..\..\01 ANSYS\07 Optimization\02_parametrisation1\02_exports\backups';
% runRef = 'hdkrkyy';   % 100-200 Hz, 1:1.14:1.9 starting point x0 = [7 8 6 4].
% runRef = 'wjiesbx';   % 100-200 Hz, 1:0.85:1.4 starting point x0 = [7 8 6 4].
% runRef = 'fecjeep';   % 100-200 Hz, 1:1.14:1.9 starting point cuboid.
% runRef = 'wxrayqc';   % 100-200 Hz, starting point ISO RR (ratios rxy = 0.9264, rxz = 1.2634).
% runRef_all = {'hdkrkyy','wjiesbx','fecjeep'};
% runRef_all = {'hdkrkyy'};
% runRef_all = {'wjiesbx'};
% runRef_all = {'fecjeep'};
runRef_all = {'snqcogt'};
runRef_all = {'iillxzc'};
runRef_all = {'tyickyz'};
% runRef_all = {'kwwcwpj'};
runRef_all = {'drmmvgu'};
% runRef_all = {'rxztxfn'};
% runRef_all = {'kwwcwpj','snqcogt','iillxzc','tyickyz','drmmvgu'};
% runRef_all = {'snqcogt','iillxzc','tyickyz','drmmvgu'};
% runRef_all = {'tyickyz','drmmvgu'};

% Robustness check to be performed
% check = 'FlowRes';
% check = 'Geometry';
check = 'SampleDims';

% BOOLEANS
dontsave = 0;       % If 1, overrides the saving of a MAT file

% Forced physical parameters (to not force, use NaN)
ds = 0.2;       % Sample thickness

rng('default')   % DEFINE RANDOM GENERATOR SEED

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

    %% PROCESS

    fname = ['RC_' runRef '_' check];
    if ~exist(['01_mat\' fname '.mat'],'file') || dontsave
        % Perform robustness checks
        status = perf_robust_checks(designData,'Check',check);    % <--- PERFORM ROBUSTNESS CHECKS

        if ~dontsave
            save(['01_mat\' fname '.mat'],'status','check','runRef')
        end
    else
        load(['01_mat\' fname '.mat'],'status','check','runRef')
    end
    
    if strcmp(check,'SampleDims')
        combinedDataForPost{idxmain}.status = status;
        combinedDataForPost{idxmain}.runRef = runRef;
    else
        %%% PLOT 1
        if isfield(designData,'workspace')
            xlims = [designData.workspace.ManualParams.fcmin/2^(1/6),...
                designData.workspace.ManualParams.fcmax*2^(1/6)-2];
        else
            xlims = [-Inf Inf];
        end

        % - - - - - - - -
        close all
        fig = figure; fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.3646 0.4861];
        hold on; grid on

        cmap = hot(length(status.RCparam));

        % Make legend
        switch check
            case 'FlowRes'
                C1    = cell(size(string(status.RCparam(:))));
                C1(:) = {' Ns/m$^4$'};
                C2    = cell(size(string(status.RCparam(:))));
                C2(:) = {'$\sigma$='};
                leg = string(C2) + string(round(status.RCparam(:))) + string(C1);
        end

        for jj = 1:size(status.alphab,2)
        % for jj = 2
            switch check
                case 'FlowRes'
                    ph(jj) = plot(exactToNormOTOBs(status.fc),status.alphab(:,jj),'o','color',cmap(jj,:));
                    fc_stairs = [exactToNormOTOBs(status.fc)/2^(1/6) exactToNormOTOBs(status.fc(end))*2^(1/6)];
                    stairs(fc_stairs, [status.alphab(:,jj); status.alphab(end,jj)],'-','color',cmap(jj,:))
                    plot(exactToNormOTOBs(status.fref),status.alpharef(:,jj),'o','color',cmap(jj,:));
                    fref_stairs = [exactToNormOTOBs(status.fref)/2^(1/6) exactToNormOTOBs(status.fref(end))*2^(1/6)];
                    stairs(fref_stairs, [status.alpharef(:,jj); status.alpharef(end,jj)],'--','color',cmap(jj,:))
                case 'Geometry'
                    if jj == 1
                        plot(exactToNormOTOBs(status.fref),status.alpharef,'ko');
                        fref_stairs = [exactToNormOTOBs(status.fref)/2^(1/6) exactToNormOTOBs(status.fref(end))*2^(1/6)];
                        stairs(fref_stairs, [status.alpharef; status.alpharef(end)],'k--')
                    end
    %             case 'SampleDims'
    %                 % Robustness of sample dimensions
    %                 if jj == 1
    %                     plot(exactToNormOTOBs(status.fref),status.alpharef,'ko');
    %                     fref_stairs = [exactToNormOTOBs(status.fref)/2^(1/6) exactToNormOTOBs(status.fref(end))*2^(1/6)];
    %                     stairs(fref_stairs, [status.alpharef; status.alpharef(end)],'k--')
    %                 end
    %                 fc_stairs = [exactToNormOTOBs(status.fc)/2^(1/6) exactToNormOTOBs(status.fc(end))*2^(1/6)];
    %                 stairs(fc_stairs, [status.alphab(:,jj); status.alphab(end,jj)],'-','color',cmap(jj,:))
            end
        end

        % TEMPORARY -- CORRECT DATA FOR BIAS INTRODUCED BY DECORRELATION
        switch runRef
            case 'wjiesbx'
                status.alphab(3,:) = status.alphab(3,:) + 0.06;
            case 'hdkrkyy'
                status.alphab(3:4,:) = status.alphab(3:4,:) + 0.04;
            case 'snqcogt'
                status.alphab(3:4,:) = status.alphab(3:4,:) + 0.04;
                status.alphab(2,:) = status.alphab(2,:) + 0.02;
            case 'kwwcwpj'
                status.alphab(1:2,:) = status.alphab(1:2,:) - 0.01;
            case 'drmmvgu'
    %             status.alphab(2,:) = status.alphab(2,:) + 0.01;
    %             status.alphab(3,:) = status.alphab(3,:) + 0.03;
    %             status.alphab(4,:) = status.alphab(4,:) + 0.03;
            case 'rxztxfn'
                status.alphab(1,:) = status.alphab(1,:) - 0.02;
                status.alphab(2,:) = status.alphab(2,:) - 0.02;
                status.alphab(3,:) = status.alphab(3,:) + 0.01;
                status.alphab(4,:) = status.alphab(4,:) + 0.02;
        end

        hold on
        switch check
            case 'Geometry'
    %             errorbar(status.fc,mean(status.alphab,2),std(status.alphab,[],2),'.','color','r');
    %             plot(status.fc,mean(status.alphab,2),'.','color','r');
    %             plot(status.fc,max(status.alphab,[],2),'.','color','r');
    %             plot(status.fc,min(status.alphab,[],2),'.','color','r');
                isout = isoutlier(status.alphab.','quartiles');
                alphab_clean = status.alphab.';
                alphab_clean(isout) = NaN;

                % Inform user
                disp(['g(v) ranges from ' num2str(min(sqrt(mean((alphab_clean - status.alpharef.').^2,2,'omitnan')))) ...
                    ' to ' num2str(max(sqrt(mean((alphab_clean - status.alpharef.').^2,2,'omitnan'))))])

                bp = boxplot(alphab_clean,'positions', exactToNormOTOBs(status.fc),...
                    'widths',[5 6 7 8],'whisker',Inf);     % BOX PLOT
                set(bp,{'linew'},{1})
        end

        switch check
            case 'FlowRes'
                legend(ph,leg,'location','se')
            case 'Geometry'
                %
        end
        OTOBticks
        xlim(xlims)
        ylim([0 Inf])

        %%% PLOT 2
        if strcmp(check,'FlowRes')

            showMin = 0;    % If 1, highlights the minimum delta on the x-axis
            deltaType = 'relative';  % Computes the relative difference det-SEA/sim
            deltaType = 'absolute';  % Computes the absolute difference det-SEA/sim

            % - - - - - - - -
            close all
            fig = figure; fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.45 0.2];
            hold on; grid on

            cmap = hot(length(status.RCparam));

            % Make legend
            switch check
                case 'FlowRes'
                    C1    = cell(size(string(status.RCparam(:))));
                    C1(:) = {' Ns/m$^4$'};
                    C2    = cell(size(string(status.RCparam(:))));
                    C2(:) = {'$\sigma$='};
                    leg = string(C2) + string(round(status.RCparam(:))) + string(C1);
            end

            delta = zeros(size(status.alphab,1),length(status.RCparam));
            for jj = 1:length(status.RCparam)
                switch deltaType
                    case 'relative'
                        delta(:,jj) = abs(status.alphab(:,jj) - status.alpharef(:,jj))./status.alpharef(:,jj)*100;
                    case 'absolute'
                        delta(:,jj) = abs(status.alphab(:,jj) - status.alpharef(:,jj));
                end
            end

            if showMin
                % Find ideal performance
                mdelta = mean(delta,1);
                [~,idxm] = min(mdelta);
                tx = sort([2000 5000 20000 30000 40000 50000 status.RCparam(idxm)]);
            else
                tx = [2000 5000 10000 20000 30000 40000 50000];
            end

            imagesc(status.RCparam,status.fc,delta)
            cb = colorbar; cb.TickLabelInterpreter = 'latex';
            switch deltaType
                case 'relative'
                    cb.Label.String = '$\Delta_\mathrm{r}\alpha$ [\%]';
                    caxis([0 20])
                case 'absolute'
                    cb.Label.String = '$\Delta_\mathrm{a}\alpha$';
                    caxis([0 0.225])
            end
            cb.Label.Interpreter = 'latex';
            axis tight; 
            xlim([(2000+4.526315789473684e+03)/2 Inf])
            shading interp
            ax = gca;
            xticks(tx)
            xticklabels(string(round(tx/1e3)))
            yticks(status.fc)
            yticklabels(string(exactToNormOTOBs(status.fc)))

            if showMin
                lab = ax.XTickLabel(ax.XTick == status.RCparam(idxm));
                ax.XTickLabel{ax.XTick == status.RCparam(idxm)} = ...
                    ['$\textbf{' lab{:} '}$'];  % Highlight minimum
            end

            ylabel '$f_\mathrm{c}$ [Hz]'
            xlabel '$\sigma$ [kNs/m$^4$]'
            ax.YScale = 'log';
            ax.TickDir = 'out';
            ax.YMinorTick = 'off';
        end
    end

    
    if 0
        if strcmp(check,'Geometry')
            fname = ['RC_' runRef '_' check];
        else
            fname = ['RC_' runRef '_' check '_' deltaType];
        end

        exportfigure(gcf, fname, '00_fig',1)
    end
    
end

%% POST PROCESS FOR SAMPLE DIMS ROBUSTNESS CHECK

if strcmp(check,'SampleDims')
    postprocess_sampledims(combinedDataForPost, ~dontsave)
end