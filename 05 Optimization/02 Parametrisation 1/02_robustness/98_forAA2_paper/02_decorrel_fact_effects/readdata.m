function data = readdata(files)
% readdata -- Reads data output after Ansys analysis for the RR design
% optimisation.
%
% >>> Inputs:
% -files [N*1 cell of string] - Files to be read.
% >>> Outputs:
% -data [struct] - Sructure data.

% (c) Paul Didier - 20-Apr-2021 11:01

% Defaults
if nargin <= 0
	%
end

% (1) Frequencies
if any(contains({files.name}, 'freqsMA.txt'))
    idx = contains({files.name}, 'freqsMA.txt');
    fid = fopen([files(idx).folder '\' files(idx).name]);
    cdata = textscan(fid, '%f %f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
    fclose(fid);
    Nfeig = size(cdata{:},1);    % Number of eigenfrequencies
    freqs = cdata{:}(:,1);    % Eigenfrequencies
else
    disp('WARNING: No eigenfrequency file detected in Ansys output (<freqsMA.txt> absent).')
end

% (2) Parameters
if any(contains({files.name}, 'paramsMA.txt'))
    idx = contains({files.name}, 'paramsMA.txt');
    fid = fopen([files(idx).folder '\' files(idx).name]);
    cdata = textscan(fid, '%f %f %f %f %f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
    fclose(fid);
    cdata = cdata{:};
    V = cdata(1,1);         % Physical parameters
    S = cdata(1,2);
    meshsizeRoom = cdata(1,3);
    c = cdata(1,4);
    rho = cdata(1,5);
    Csx = cdata(2,1);           % Sample parameters
    Csy = cdata(2,2);           % Sample parameters
    Lsx = cdata(2,3);           % Sample parameters
    Lsy = cdata(2,4);           % Sample parameters
    ds = cdata(2,5);            % Sample parameters
    s_coords = cdata(3:end,1:3);        % Sources coordinates
    % Create substructure
    params_MA = struct('V',V,'S',S,'e_room',meshsizeRoom,'c',c,'rho',rho,...
        's_coords',s_coords,'Csx',Csx,'Csy',Csy,'Lsx',Lsx,'Lsy',Lsy,'ds',ds);
else
    disp('WARNING: No Modal Analysis parameters file detected in Ansys output (<paramsMA.txt> absent).')
end

% (3) Source mode shapes
if any(contains({files.name}, 'psMA'))
    idx = contains({files.name}, 'psMA');
    tmp = 1:length(files);
    idx = tmp(idx);
    for ii = 1:length(idx)   % Adapted to multiple-sources export -- 20210428
        fid = fopen([files(idx(ii)).folder '\' files(idx(ii)).name]);
        cdata = textscan(fid, '%f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
        fclose(fid);
        ps(:,ii) = cdata{:};  % Source mode shapes
    end
else
    disp('WARNING: No source mode shapes file detected in Ansys output (<psMA_XX.txt> absent).')
end

% (4) Mode shapes on sample surface
% Read coordinates
if any(contains({files.name}, 'sampleMSc.csv'))
    idx = contains({files.name}, 'sampleMSc.csv');
    data = csvread([files(idx).folder '\' files(idx).name]);
    nNodes = length(data)/3;
    coords = [data(1:nNodes), data(nNodes+1:2*nNodes), data(2*nNodes+1:end)];
else
    disp('WARNING: No sample surface coordinates file detected in Ansys output (<sampleMSc.csv> absent).')
end
% Read mode shapes
if any(contains({files.name}, 'sampleMSv.csv'))
    idx = contains({files.name}, 'sampleMSv.csv');
    data = csvread([files(idx).folder '\' files(idx).name]);
    nModes = length(data)/nNodes;
    m = zeros(nNodes,4,nModes);
    for jj = 1:nModes
        m(:,4,jj) = data(jj:nModes:end);
        m(:,1:3,jj) = coords;
    end
else
    disp('WARNING: No sample surface mode shapes file detected in Ansys output (<sampleMSv.csv> non-existent).')
end

data = struct('params',params_MA,'f',freqs,'ps',ps,'nModes',nModes,'MS',m); 
  
end