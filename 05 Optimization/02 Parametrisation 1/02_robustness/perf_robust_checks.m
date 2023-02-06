function status = perf_robust_checks(designData,varargin)
% perf_robust_checks -- Wrapping function to perform robustness checks
% (RCs) given a certain RR design's FEM-simulated data.
%
% >>> Inputs:
% -designData [struct] - Design data read and arranged via <readdesigndata()>. 
% -varargin:
%   -'Check' [string] - Check to be performed.
%       -If 'FlowRes': RC across sample flow resistivities (default).
%       -If 'Geometry': RC across slight geometrical changes.
% >>> Outputs:
% -status [struct] - Outcome of RC check.

% (c) Paul Didier - 11-May-2021 10:21


% Check varargin
check = 'FlowRes';
idx_varargin = 1:length(varargin);
if any(strcmp(varargin, 'Check'))
    check = varargin{idx_varargin(strcmp(varargin, 'Check')) + 1};
end

% Defaults
switch check
    case 'FlowRes'
        delete('uncorrel_factor.mat');
        status = perf_flowres_RC(designData);
    case 'Geometry'
        status = perf_geom_RC(designData);
end

end