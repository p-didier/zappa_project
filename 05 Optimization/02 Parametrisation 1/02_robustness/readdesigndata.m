function data = readdesigndata(designDataPath,varargin)
% readdesigndata -- Reads data in a folder containing the data necessary to
% assess a certain RR design.
%
% >>> Inputs:
% -designDataPath [string] - Path to folder containing design data.
% -varargin:
%   -'Workspace' [bool] - If true, fetches the workspace too (default:
%   false).
%   -'OptiOutcome' [bool] - If true, fetches the final optimization run
%   outcome too (default: false).
% >>> Outputs:
% -data [struct] - Data read and arranged. 

% (c) Paul Didier - 11-May-2021 10:20

% Interpret varargin
optioutflag = 0; wsflag = 0;
idx_varargin = 1:length(varargin);
if any(strcmp(varargin, 'OptiOutcome'))
    optioutflag = varargin{idx_varargin(strcmp(varargin, 'OptiOutcome')) + 1};
end
if any(strcmp(varargin, 'Workspace'))
    wsflag = varargin{idx_varargin(strcmp(varargin, 'Workspace')) + 1};
end

% -------- Get Ansys outputs and Matlab post-processing data for iteration
currdir = pwd;
cd(designDataPath)

% Read Modal Analysis data
% (1) Frequencies
fid = fopen('freqsMA.txt');
cdata = textscan(fid, '%f %f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
fclose(fid);
Nfeig = size(cdata{:},1);    % Number of eigenfrequencies
freqs = cdata{:}(:,1);    % Eigenfrequencies

% (2) Parameters
fid = fopen('paramsMA.txt');
cdata = textscan(fid, '%f %f %f %f %f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
fclose(fid);
cdata = cdata{:};
Vroom = cdata(1,1);         % Physical parameters
Sroom = cdata(1,2);         % Physical parameters
meshsizeRoom = cdata(1,3);  % Physical parameters
c = cdata(1,4);             % Physical parameters
rho = cdata(1,5);           % Physical parameters
Csx = cdata(2,1);           % Sample parameters
Csy = cdata(2,2);           % Sample parameters
Lsx = cdata(2,3);           % Sample parameters
Lsy = cdata(2,4);           % Sample parameters
ds = cdata(2,5);            % Sample parameters
s_coords = cdata(3:end,1:3);    % Source coordinates
% Create substructure
params_MA = struct('V',Vroom,'S',Sroom,'e_room',meshsizeRoom,'c',c,'rho',rho,...
    'Csx',Csx,'Csy',Csy,'Lsx',Lsx,'Lsy',Lsy,'ds',ds,...
    's_coords',s_coords);

% (3) Source pressure
files = dir; idx = [];
for ii = 1:length(files)
    if length(files(ii).name) > 2 && any(contains(files(ii).name,'psMA'))
        idx(end+1) = ii;
    end
end
for ii = 1:length(idx)  
    fid = fopen(files(idx(ii)).name);
    cdata = textscan(fid, '%f', 'delimiter', ',', 'HeaderLines', 0, 'CollectOutput',1);
    fclose(fid);
    ps(:,ii) = cdata{:};  % Real pressure
end

% (4) Post-processing data
matfile = files(contains({files.name},'.mat'));
postProcData = load(matfile.name);

% (5) Mode shapes on sample surface
% Read coordinates
if exist('sampleMSc.csv', 'file')
    data = csvread('sampleMSc.csv');
else
    error('No data file or incorrect name for mode shape nodal coordinates.')
end
nNodes = length(data)/3;
coords = [data(1:nNodes), data(nNodes+1:2*nNodes), data(2*nNodes+1:end)];
% Read mode shapes
if exist('sampleMSv.csv', 'file')
    data = csvread('sampleMSv.csv');
else
    error('No data file or incorrect name for mode shapes.')
end
nModes = length(data)/nNodes;
m = zeros(nNodes,4,nModes);
for jj = 1:nModes
    m(:,4,jj) = data(jj:nModes:end);
    m(:,1:3,jj) = coords;
end

% Extract just folder name for referencing
idxslash = regexp(designDataPath,'\');
foldername = designDataPath(idxslash(end)+1:end);

data = struct('params',params_MA,'f',freqs,'ps',ps,...
    'nModes',nModes,'MS',m,...
    'postProcData',postProcData,...
    'fullpath',designDataPath,'foldername', foldername); 

cd(currdir)

% Save run reference label
runref = designDataPath(end-10:end-4);
data.('runref') = runref;

if optioutflag
    % -------- Get final optimization run outcome
    % Fetch relevant data
    username = getenv('USERNAME');
    pathToMat = ['C:\Users\' username '\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\05 Optimization\02 Parametrisation 1\01_mat\02_outcomes'];
    files = dir(pathToMat); 
%     if size(files,1) == 0
%         pathToMat = 'C:\Users\paulr\Dropbox\BELGIUM\KUL\Research\04 Simulations\05 Optimization\02 Parametrisation 1\01_mat\02_outcomes';
%         files = dir(pathToMat); 
%     end
    files(1:2) = [];
    files = files(contains({files.name},runref));
    if size(files,1) == 0
        data.('opti_outcome') = 'NO FILE FOUND';
    else
        load([files.folder '\' files.name],'res')
        % Export it
        data.('opti_outcome') = res;
    end
end

if wsflag
    % -------- Get final workspace data
    username = getenv('USERNAME');
    pathToMat = ['C:\Users\' username '\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\05 Optimization\02 Parametrisation 1\01_mat\01_workspaces'];
    files = dir(pathToMat); 
%     if size(files,1) == 0
%         pathToMat = 'C:\Users\paulr\Dropbox\_BELGIUM\KUL\Research\04 Simulations\05 Optimization\02 Parametrisation 1\01_mat\01_workspaces';
%         files = dir(pathToMat); 
%     end
    files(1:2) = [];
    files = files(contains({files.name},runref));
    if size(files,1) == 0
        data.('workspace') = 'NO FILE FOUND';
    else
        load([files.folder '\' files.name],'mp','options','rd','x0')
        ws = struct('ManualParams',mp,'OptimOptions',options,'RoomDim',rd,'InitGuess',x0);
        data.('workspace') = ws;
    end
end

end