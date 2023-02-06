function dims = find_valid_sample_dimensions(n)
% find_valid_sample_dimensions -- Finds a number of ISO 354-valid sample
% dimensions (surface area between 10 and 12 squared meters, and length-
% to-width ratio between 0.7 and 0.9 [ISO/CD 354:2019]).
%
% >>> Inputs:
% -n : int
%    Number of possible XY dimensions combinations to find.
% >>> Outputs:
% -dims : [n x 2] array (float)
%    List of dimensions [m].

% (c) Paul Didier - 06-Feb-2023 12:35
% SOUNDS ETN - KU Leuven ESAT STADIUS
% ------------------------------------

% Defaults
if nargin <= 0
	%
end

rng('default')

% ISO 354 guidelines
surfAreaBounds = [10, 12];  % [m^2]
ltwRatioBounds = [0.7, 0.9];

dims = zeros(n, 2);
for ii = 1:n
    % Set a random surface area
    surfaceArea = (max(surfAreaBounds) - min(surfAreaBounds)) * ...
        rand(1, 1) + min(surfAreaBounds);
    % Set a random length-to-width ratio
    ltwRatio = (max(ltwRatioBounds) - min(ltwRatioBounds)) * ...
        rand(1, 1) + min(ltwRatioBounds);
    % Derive corresponding dimensions
    w = sqrt(surfaceArea / ltwRatio); 
    l = w * ltwRatio;
    dims(ii, :) = [l, w];
end

end