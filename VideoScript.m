%% Simulation Video 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: VideoScript.m
%
% Description: create a video illustrating simulation results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Video

hFig = figure;
hold on

% Video object options
vobj = VideoWriter('movie.mp4', 'MPEG-4');
vobj.FrameRate = 16;
vobj.Quality = 100;
open(vobj);

% Resample the time vector in order to obtain equally spaced samples and 
% eliminate duplicates due to jumps
idxInit = 1;
tk = t(idxInit); 
steps = 2000;
deltaT = 80*options.MaxStep;
if deltaT < options.MaxStep
   error('VideoScript ERROR: deltaT must be larger than simulator MaxStep')
end
k = 1;
timeVec = 0;
for idx = idxInit:length(t)
    if abs(t(idx) - tk) <= deltaT
        timeVec(k) = idx;
        k = k+1;
        tk = tk + deltaT;
    end
end

% Scaling factor for text and graphic objects
scale = 1.5;
%%%%%%
% Note: always include 'hold on' and 'hold off' in each subplot
%%%%%%
Ti = 1;
for idx = timeVec
    cla(gca)
    % Script to plot single frame
    PlotFrame;
    % Figure properties
        % Axes font size
    set(findall(hFig,'type','axes'),'fontsize',scale*12)
        % Text font size
    set(findall(hFig,'type','text'),'fontSize',scale*14)
        % Figure dimension
    set(hFig,'units','pixels','position',[100, 0, 860, 600])
        % Background color
    set(hFig,'color','w')
    % Get current frame
    F = getframe(hFig);
    % Write frame to video object
    writeVideo(vobj, F);  
end

% Close video object
close(vobj);