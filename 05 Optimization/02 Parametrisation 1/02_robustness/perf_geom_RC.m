function status = perf_geom_RC(dd)
% perf_geom_RC -- Performs a robustness check (RC) on a given RR design
% across slight building tolerances and neighbouring geometries.
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
inpath = 'C:\Users\u0137935\Dropbox\BELGIUM\KU Leuven\Research\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\01_input_files';
outpath = 'C:\Users\u0137935\Dropbox\BELGIUM\KU Leuven\Research\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1\02_exports';
   
% -------------------------
tolrange = 0.05;    % Tolerance range on coordinates [m]

% Identify variables that were used to build the room
niter = str2double(dd.foldername(end-2:end));
if ~ischar(dd.opti_outcome)
    itervars = dd.opti_outcome.varevol(:,niter);
else
    if strcmp(dd.runref,'hdkrkyy')    % Case-specific -- 20210513
        itervars = [7 8 4.06 5.5];
    end
end

N = 50;

% % % Generate RC coordinate sets
% % coords = getWiggledCoords(itervars,tolrange);
% % % Select N random sets
% % coords = coords(randi(size(coords,1),[N,1]),:);

% Generate RC coordinate sets
coords = getTolCoords(itervars,tolrange,N);

% % % TEMPORARY
% % coords = [7,8.05000000000000,4.11000000000000,5.55000000000000];

% Get relevant workspace data
rxy = dd.workspace.ManualParams.rxy;
rxz = dd.workspace.ManualParams.rxz;
V = dd.workspace.ManualParams.V;
fcmin = dd.workspace.ManualParams.fcmin;
fcmax = dd.workspace.ManualParams.fcmax;
mof = dd.workspace.ManualParams.mof;

% Interpret coordinates
for ii = 1:size(coords,1)
    tic
    %%% Step 1 - Write the input file for Ansys
    params = interpret_coords(coords(ii,:),rxy,rxz,V);
    params.('fmin') = fcmin/2^(1/6)/mof;
    params.('fmax') = fcmax*2^(1/6)*mof;

    % TMP -- QUICK TESTS
% % %     params.('fmin') = 20; params.('fmax') = 50; dd.workspace.ManualParams.flowres_ds= [28e3,0.2];
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

%% 

function coords = getWiggledCoords(vars,tolrange)
% Generates all possible combinations of wiggles around the central values
% <vars> by a distance <tolrange> -- e.g. for a single variable <x>, the
% function <getWiggledCoords> returns <[x - tolrange; x + tolrange]>.

wiggles = cell(length(vars),1);
for ii = 1:length(vars)
    wiggles{ii} = [vars(ii) - tolrange; vars(ii); vars(ii) + tolrange];
end
coords = allcomb(wiggles{:});

coords((size(coords,1) - 1)/2 + 1,:) = [];  % Eliminate the non-wiggled combination.

stopline = 1;

end

function coords = getTolCoords(vars,tolrange,N)
% Generates <N> randomly shifted sets of coordinates <vars> within a certain
% range <tolrange>.

coords = zeros(N, length(vars));

for ii = 1:N
    coords(ii,:) = vars + (2*rand(size(vars)) - 1)*tolrange;
end

stopline = 1;

end