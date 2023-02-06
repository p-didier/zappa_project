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
   
% ------------GLOBAL VARIABLES-------------
nDimCombs = 3;
% -----------------------------------------

% Set sample dimensions to test
sampleXYdims = find_valid_sample_dimensions(nDimCombs);
if ~any(sampleXYdims == [3.0, 3.6])
    % Make sure the original (3 x 3.6) sample dimensions are included
    sampleXYdims = [[3.0, 3.6]; sampleXYdims(1:end-1, :)];
end

% Identify variables that were used to build the room
niter = str2double(dd.foldername(end-2:end));
itervars = dd.opti_outcome.varevol(:,niter);

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

% Loop over sample dimensions
for ii = 1:nDimCombs
    tic
    %%% Step 1 - Write the input file for Ansys
    params.('Lsx') = sampleXYdims(ii, 1);
    params.('Lsy') = sampleXYdims(ii, 2);
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
    [~,results] = getSNV(FEdata,dd.workspace.ManualParams,0);
    
    % Compile results for export
    alphab_all(:,ii) = results.alphasimb;
    toc
    stopline = 1;
end

status = struct('alphab',alphab_all,'alpharef',results.alpharef,...
    'fc',results.fc,'fref',results.fref,'RCparam',coords);

end