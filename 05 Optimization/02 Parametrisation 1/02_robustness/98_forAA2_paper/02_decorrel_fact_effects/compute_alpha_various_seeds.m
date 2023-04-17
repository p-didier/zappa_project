function results = compute_alpha_various_seeds(FEdata, mp, nRuns, testType)
% compute_alpha_various_seeds -- Compute absorption coefficient given Ansys
% outputs. Using the FE-RR way to estimate ACs. Computing the modal force
% vector using various seeds for the multi-source decorrelation factor.
%
% >>> Inputs:
% -FEdata [struct] - Simulation output organised in a structure.
% -mp [struct] - Manually-input parameters.
% -nRuns [int] - Number of (Monte-Carlo) runs.
% -testType [str] - 'changing_seed' or 'changing_state'
% >>> Outputs:
% -SNV [float, -] - Single-Number-Value associated to input parameters.
% -plotdata [struct] - Data useful for subsequent dynamic plotting.

% (c) Paul Didier - 20-Apr-2021 10:54

%% INTERPRET FEM DATA

% Extract data from FEM exports
sd(1) = FEdata.params.Lsx;      % Sample length.
sd(2) = FEdata.params.Lsy;      % Sample width.
c = FEdata.params.c;            % Speed of sound in air.
rho = FEdata.params.rho;        % Air density.
V = FEdata.params.V;            % Air volume. 

%% INTERPRET MANUALLY INPUT PARAMETERS

Xi = mp.flowres_ds(1);          % Effective flow resistivity, used to get surface impedance <Zs>.
ds_eff = mp.flowres_ds(2);      % Effective thickness, used to get surface impedance <Zs>.
fcmin = mp.fcmin;               % Low bound on center-frequencies to consider.
fcmax = mp.fcmax;               % High bound on center-frequencies to consider.
Tempty = mp.Tempty;             % Empty RT.
dg = mp.dg;                     % 2-D numerical integration spatial grid resolution.

fmin = fcmin/2^(1/6)/mp.mof;
fmax = fcmax*2^(1/6)*mp.mof;

%% PRE-PROCESSING

% Compute coupling matrix
disp('Computing coupling matrix...')
Cpp = getCouplingMat(FEdata.MS,Inf,dg,'silent',1);
% Extract eigenfrequencies
feig = FEdata.f; 
idxeig = findIndex(feig, fmin);
% Construct the frequency vector,
freq = getFreqVector(fmin,fmax,mp.nfpb);
% ... and the room's damping loss factor
eta = getEtaFromTempty(Tempty,freq);

%% OBTAIN SIMULATED ABSORPTION COEFFICIENT AND REFERENCE DET-SEA VALUE

% Compute surface admittance
Zs = Z_Miki(rho,c,freq,Xi,ds_eff,0); 
beta = rho*c./Zs;

if strcmp(testType, 'changing_state')
    rng('default')
else
    seeds = 1:nRuns;  % get the seeds
end

for ii = 1:length(nRuns)

    disp(['Running with rng(' num2str(seeds(ii)) ')... [run ' num2str(ii) '/' num2str(length(seeds)) ']'])
    
    if strcmp(testType, 'changing_seed')
        % Set seed
        rng(seeds(ii))
    end

    % Compute the modal force (force in spatial domain, uncorrelated sources)
    fvec = getModalForce(FEdata.ps(idxeig:end,:),freq,true);   
    
    % Compute modal coordinatess
    disp('Computing modal coordinates...')
    [qsample,qempty] = getModalCoords(feig(idxeig:end),freq,...
                            Cpp(idxeig:end,idxeig:end),fvec,beta,eta,c);  
    
    % Get absorption coefficients from energy method
    disp('Computing absorption coefficients...')
    [~,~,alpha] = getAC_fromEnergy(qempty,qsample,freq,feig(idxeig:end),fvec.',...
                                    V,prod(sd),c,rho);
    
    % Integrate results over bands
    [alphab,fc,~,fu] = harmToBands(alpha,freq,3);
    % Discard too low bands
    idx = fc < fcmin/2^(1/6);
    alphab(idx) = []; fc(idx) = []; fu(idx) = [];
    % Discard too high bands
    idx = fu > fcmax*2^(1/6);
    alphab(idx) = []; fc(idx) = [];
    
    % Fetch reference data
    [alpharef,fref] = getOFCTtarget(Xi,ds_eff,sd);
    
    % Find overlapping bands between test and reference
    idxovlapref = false(length(fc),1);
    idxdiscard = false(length(fc),1);
    for jj = 1:length(fc)
        if any(round(exactToNormOTOBs(fref)) == round(exactToNormOTOBs(fc(jj))))
            disp(['Target value found for fc = ' num2str(exactToNormOTOBs(fc(jj))) ' Hz. Using result.'])
            idxovlapref(findIndex(fref,fc(jj))) = true;
        else
            disp(['No target value for fc = ' num2str(exactToNormOTOBs(fc(jj))) ' Hz. Discarding result.'])
            idxdiscard(jj) = true;
        end
    end
    alphab(idxdiscard) = []; % No reference value for these bands
    fc(idxdiscard) = []; % No reference value for these bands
    alpharef = alpharef(idxovlapref);
    fref = fref(idxovlapref);

    alpha_all(:, ii) = alpha;
    alphab_all(:, ii) = alphab;
    alpharef_all(:, ii) = alpharef;
end
    
% Prepare data export for subsequent dynamic plotting
results.fc = fc;
results.f = freq;
results.alpha_all = alpha_all;
results.alphab_all = alphab_all;
results.alpharef_all = alpharef_all;
results.fref = fref;
results.Ns = size(FEdata.ps,2);
results.rxy = mp.rxy;
results.rxz = mp.rxz;
results.V = mp.V;

end