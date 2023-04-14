clear; close all; clc;

% -- Purpose of script
% 

% (c) Paul Didier - 23-Mar-2023 09:54
% SOUNDS ETN - KU Leuven ESAT STADIUS

%% INIT

rng('default');
count = 1;
while 1
    state = rng;
    plot(state.State)
    title(count)
    hold off
%     state.State(end)
%     pause(1)
    a = rand(1);
    count = count + 1;
    if state.State(end) == 2
        disp(state.State(end - 1))
    end
end

%% PROCESS


%% PLOT

% close all
% fig = figure; fig.Units = "Normalized"; fig.Position = [0.3177 0.3958 0.3646 0.4861];
% hold on; grid on

% ax = gca;

%% FUNCTIONS

