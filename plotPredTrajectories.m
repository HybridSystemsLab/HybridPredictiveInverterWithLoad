%% Plot predicted Trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: plotPredTrajectories.m
%
% Description: plot of all predicted trajectories at each jump.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Cap omega H epsTilt rhoStar deltaBar
colors = lines;
cla;
% Plot reference zR
[~, xR] = ode45(@(t, xR) f_reference(xR), TPspan, [iR, vR]);

% Plot solutions
for i = 1:length(q0)
    if q0(i) == qplus
        xP = solPlot{i}.x;
        P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
        for k = 1:size(xP,1)
            patchell(P/deltaBar, [xP(k,4), xP(k,5)],'patch',colors(6,:));
        end
        for k = 1:size(xP,1)
            patchell(P/rhoStar, [xP(k,4), xP(k,5)],'patch',colors(5,:));
        end
    end
end

colSol = [colors(1,:); colors(3,:); colors(4,:)];
for i = 1:length(q0)
    h(i) = plot(solPlot{i}.x(1:end,2), solPlot{i}.x(1:end,3), ...
        'color', colSol(i,:), 'LineWidth', 2);
end

% Plot markers
for i = 1:length(q0)
    xqplus = get(h(i), 'XData');
    yqplus = get(h(i), 'YData');
    if i == find(q0 == qplus)
        % Star marker (qplus)
        hqplus = scatter(xqplus(1,len(i)), yqplus(1,len(i)), 300 ,'r', ...
            'filled', 'p', 'MarkerEdgeColor', 'k');
    else
        % Diamond markers
        scatter(xqplus(1,len(i)), yqplus(1,len(i)), 150 , 'k', ...
            'filled', 'd', 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'k');
    end
end

legend([h,hqplus], {'q = -1', 'q = 0', 'q = 1','ball', 'qplus'})
set(gca, 'FontSize',12)
