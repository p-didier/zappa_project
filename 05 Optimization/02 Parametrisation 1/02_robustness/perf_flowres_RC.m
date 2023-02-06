function status = perf_flowres_RC(dd)
% perf_flowres_RC -- Performs a robustness check (RC) on a given RR design
% across a range of sample flow resitivities.
%
% >>> Inputs:
% -dd [struct] - Design data read and arranged via <readdesigndata()>. 
% >>> Outputs:
% -status [] - .

% (c) Paul Didier - 11-May-2021 10:27

% Hard-coded
Tempty = 10;  % Empty room RT [s]

% ------

% Flow resistivities to test
flowResRange = linspace(2,50,20)*1e3;
flowResRange(1) = [];

% Interpret inputs
sd = [dd.params.Lsx, dd.params.Lsy];
ds = dd.params.ds;
Cpp = dd.postProcData.Cpp;
ap = dd.postProcData.plotdata;
freq = ap.f;
feig = dd.f; 
c = dd.params.c;
rho = dd.params.rho;
V = dd.params.V;

% Damping loss factor
eta = getEtaFromTempty(Tempty,freq);
% Modal force (force in spatial domain, uncorrelated sources)
idxeig = findIndex(feig, min(freq));
fvec = getModalForce(dd.ps(idxeig:end,:),freq);     

% For each flow resisitivity...
for ii = 1:length(flowResRange)
    disp(['Computing RC for flowres #' num2str(ii) '\' num2str(length(flowResRange)) '...'])
    flowRes = flowResRange(ii);
    
    % Fetch reference data
    [alpharef,fref] = getOFCTtarget(flowRes,ds,sd,'RC',1);
    
    % --------- Compute measured absorption coefficient
    Zs = Z_Miki(rho,c,freq,flowRes,dd.params.ds,0);   
    beta = rho*c./Zs;  % Surface admittance

    % Compute modal coordinates
    [qsample,qempty] = getModalCoords(feig(idxeig:end),freq,...
                            Cpp(idxeig:end,idxeig:end),fvec,beta,eta,c);  
    % Get absorption coefficients from energy method
    [~,~,alpha] = getAC_fromEnergy(qempty,qsample,freq,feig(idxeig:end),fvec.',...
                                    V,prod(sd),c,rho);
    
    % Moving-average smoothing to counteract the random decorrelation
    alpha = movmean(alpha,20);         % <--- SPECIFIC TO THIS SCRIPT!!
                                
    % Integrate results over bands
    [alphab,fc,~,fu] = harmToBands(alpha,freq,3);
    % Discard too low bands
    idx = fc < dd.workspace.ManualParams.fcmin/2^(1/6); alphab(idx) = []; fc(idx) = []; fu(idx) = [];
    % Discard too high bands
    idx = fu > dd.workspace.ManualParams.fcmax*2^(1/6); alphab(idx) = []; fc(idx) = [];
    

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
    
    if ~all(size(alphab) == size(alpharef))
        alpharef = alpharef.';
    end
    
    % Compile results outside of loop
    alphab_all(:,ii) = alphab;
    alpharef_all(:,ii) = alpharef;
    
    % Inform user
    rdist = abs(alphab - alpharef)./alpharef*100;
    disp(['For Xi = ' num2str(flowRes) ' Ns/m^4, rel. distance to target [%] per bands from '...
        num2str(fref(1)) ' to ' num2str(fref(end)) ' Hz:'])
    disp(rdist)
end

status = struct('alphab',alphab_all,'alpharef',alpharef_all,'fc',fc,...
    'fref',fref,'RCparam',flowResRange);

end