clear; close all; clc;
addpath(genpath('..\..\..\02 MATLAB\09 Reverberation Time sims\10 RTana_v3\toolbox'))
addpath(genpath('..\..\..\02 MATLAB\basefunctions'))

% -- Purpose of script
% Conduct robustness check over different absorbers for a cuboidal room.

% (c) Paul Didier - 28-Jun-2021 10:36
% ZAPPA

%% INIT

V = 250;
rxy = 1.14; rxz = 1.9;
rxy = 0.85; rxz = 1.4;

flowResRange = linspace(2,50,20)*1e3;
flowResRange(1) = [];

fmin = 100/2^(1/6);
fmax = 200*2^(1/6);
sd = [3 3.6]; 
ds = 0.2; 
Tempty = 10;
nxyz_max = 20;
nfpb = 1572;

%% PROCESS

res = zeros(length(flowResRange),length(noctfr(3,fmin,fmax/2^(1/6),'exact')));
for ii = 1:length(flowResRange)
    Xi = flowResRange(ii);

%     [alpharef,fref] = getOFCTtarget(Xi,ds,sd,'RC',1);
    [alpharef,fref] = getOFCTtarget(Xi,ds,sd);
    
    disp(['Computing for Xi = ' num2str(Xi/1e3) ' kNs/m4 (' num2str(ii)...
        '/' num2str(length(flowResRange)) ')...'])

    [~,~,~,alphaBands,fc] = getOFCTval(rxy,rxz,V,fmin,fmax,sd,Xi,ds,Tempty,nxyz_max,alpharef,fref,nfpb,...
        'sposType','verycorners','sDecoLoad',1);
    
    alpharef_relevant = zeros(length(fc),1);
    for jj = 1:length(fc)
        alpharef_relevant(jj) = alpharef(exactToNormOTOBs(fref) == exactToNormOTOBs(fc(jj)));
    end
    
    res(ii,:) = abs(alphaBands - alpharef_relevant);
end
disp('Done.')

%% PLOT

close all
fig = figure; fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.45 0.2];
hold on; grid on
imagesc(flowResRange,fc,res.')
cb = colorbar; cb.TickLabelInterpreter = 'latex';
cb.Label.String = '$\Delta_\mathrm{a}\alpha$';
caxis([0 0.225])
cb.Label.Interpreter = 'latex';
axis tight; 
xlim([(2000+4.526315789473684e+03)/2 Inf])
shading interp
ax = gca;
tx = [2000 5000 10000 20000 30000 40000 50000];
xticks(tx)
xticklabels(string(round(tx/1e3)))
yticks(fc)
yticklabels(string(exactToNormOTOBs(fc)))
ylabel '$f_\mathrm{c}$ [Hz]'
xlabel '$\sigma$ [kNs/m$^4$]'
ax.YScale = 'log';
ax.TickDir = 'out';
ax.YMinorTick = 'off';

%% EXPORT

rxy = 0.85; rxz = 1.4;

if 0
    exportfigure(gcf, ['out_V' num2str(V) 'rxy' num2str(rxy*100)...
        'rxz' num2str(rxz*100)], '00_fig\02_cuboids')
end

