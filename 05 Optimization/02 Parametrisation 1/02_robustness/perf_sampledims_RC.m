function status = perf_sampledims_RC(dd)
% perf_sampledims_RC -- Performs a robustness check (RC) on a given RR 
% design across different sample XY (2D) dimensions.
%
% >>> Inputs:
% -dd [struct] - Design data read and arranged via <readdesigndata()>. 
% >>> Outputs:
% -status [] - .

% (c) Paul Didier - 12-May-2021 09:28

% Defaults
if nargin <= 0
	%
end

% Paths
% inpath = 'C:\Users\u0137935\Dropbox\BELGIUM\KU Leuven\Research\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\01_input_files';
% outpath = 'C:\Users\u0137935\Dropbox\BELGIUM\KU Leuven\Research\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\02_exports';
inpath = 'C:\Users\u0137935\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\01_input_files';
outpath = 'C:\Users\u0137935\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\02_exports';

% Set sample dimensions to test
load sampleDims_robCheck
sampleXYdims = sampleDims;
if ~any(sampleXYdims == [3.0, 3.6])
    % Make sure the original (3 x 3.6) sample dimensions are included
    sampleXYdims(12, :) = [];  % get rid of problematic sample (see SOUNDS PhD journal 2023 week 09 THU)
    sampleXYdims = [[3.0, 3.6]; sampleXYdims];
end
sampleXYdims = [[3.0, 3.6]];   % DEBUGGING 20.03.2023

%%% TMP: FOR DEBUGGING    -- 23.02.2023 14:49
% sampleXYdims = sampleXYdims(1:2, :);
%%% TMP: FOR DEBUGGING    -- 02.03.2023 10:41
% sampleXYdims = [[2.90039576600427,4.10304146667245]];

% Identify variables that were used to build the room
niter = str2double(dd.foldername(end-2:end));
itervars = dd.opti_outcome.varevol(:,niter);

% Set the correct random number generator state
rngState = dd.opti_outcome.results.output.rngstate; % state of random number generator
nRandCalls = rngState.State(end) / 2;  % See AA2 paper notes, pp.11-12.
offsetRandCalls = nRandCalls - niter +69;
% Set correct seed
rng(rngState.Seed)
% Run rand() enough times to reach the correct state
for ii = 1:offsetRandCalls
    rand(1);
end

% Get relevant workspace data
rxy = dd.workspace.ManualParams.rxy;
rxz = dd.workspace.ManualParams.rxz;
V = dd.workspace.ManualParams.V;
fcmin = dd.workspace.ManualParams.fcmin;
fcmax = dd.workspace.ManualParams.fcmax;
mof = dd.workspace.ManualParams.mof;


params = interpret_coords(itervars,rxy,rxz,V);
params.('fmin') = fcmin/2^(1/6)/mof;
params.('fmax') = fcmax*2^(1/6)*mof;

nDimCombs = size(sampleXYdims, 1);
% Loop over sample dimensions
for ii = 1:nDimCombs
    tic
    %%% Step 1 - Write the input file for Ansys
    params.('Lsx') = sampleXYdims(ii, 1);
    params.('Lsy') = sampleXYdims(ii, 2);
    disp(['>>>>>> Sample dims ' num2str(ii) '/' num2str(nDimCombs) ...
        ': ' num2str(round(params.Lsx * 100)) 'x' ...
        num2str(round(params.Lsy * 100))])
    write_input_file(params,inpath,0)
    %%% Step 2 - Launch Ansys analysis
    launchflag = true;
    while launchflag
        launch_analysis(inpath,outpath)
        % Dynamically check for Ansys output
        status = detect_FE_data(outpath);
        if ~any(strcmp(status, {'LicenseError', 'FatalError'}))
            launchflag = false;
        else
            fprintf('\nError detected in Ansys prompt -- Relaunching analysis.\n')
        end
    end
    %%% Step 3 - Read and process data
    FEdata = read_FE_data(outpath,@readdata,0);
    %%% Step 4 - Post-process to derive ACs
    [~,results] = getSNV(FEdata,dd.workspace.ManualParams,0,'','',NaN);
    
    % Compile results for export
    alphab_all(:,ii) = results.alphasimb;
    alpharef_all(:,ii) = results.alpharef;
    toc
end

status = struct('alphab',alphab_all,'alpharef',alpharef_all,...
    'fc',results.fc,'fref',results.fref,'RCparam',sampleXYdims);

end