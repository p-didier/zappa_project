function status = run_with_various_seeds(dd, seeds)
% run_with_various_seeds -- Runs a full ANSYS analysis, then post-processes
% its results using different random number generator seeds, to
% investigate the impact of the randomly generated source decorrelation
% factor on the absorption coefficients estimates.
%
% >>> Inputs:
% -dd [struct] - Design data read and arranged via <readdesigndata()>.
% -seeds [list] - Seeds to consider.
% >>> Outputs:
% -Ø

% (c) Paul Didier - 23-Mar-2023 11:33
% SOUNDS ETN - KU Leuven ESAT STADIUS
% ------------------------------------

% Paths
inpath = 'C:\Users\u0137935\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\01_input_files';
outpath = 'C:\Users\u0137935\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\02_exports';

% Set sample dimensions to test
sampleXYdims = [3.0, 3.6];

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

%%% Step 1 - Write the input file for Ansys
params.('Lsx') = sampleXYdims(1);
params.('Lsy') = sampleXYdims(2);
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
%%% Step 3 - Read data
FEdata = read_FE_data(outpath,@readdata,0);
%%% Step 4 - Post-process to derive ACs
results = compute_alpha_various_seeds(FEdata,dd.workspace.ManualParams,seeds);

status = struct('alphab',results.alphab_all,'alpharef',results.alpharef_all,...
    'fc',results.fc,'fref',results.fref,'RCparam',sampleXYdims);

end