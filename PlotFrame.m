%% Plot Video Frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: PlotFrame.m
%
% Description: plot of a video frame given a time vector index idx.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Cap omega H  rhoStar epsTilt deltaBar GammaLim

xLim = 1.1*GammaLim*Cap*omega;
yLim = 1.1*GammaLim;

%% Phase plane
subplot(3,2,[1,3,5])
cla(gca)
hold on
colors = lines;
% Estimated region of attraction 
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
patchell(P/deltaBar, [iR(idx),vR(idx)], 'patch', colors(5,:));
% State z = (iL, vC)
plot(iL(Ti:idx), vC(Ti:idx), 'color', colors(1,:),'LineWidth', scale*.9);
scatter(iL(idx), vC(idx), scale*50, colors(1,:), 'filled');
% Reference zB = (iB, vB)
plot(iR(Ti:idx), vR(Ti:idx), 'color', colors(2,:),'LineWidth', scale);
scatter(iR(idx), vR(idx), scale*50, colors(2,:), 'filled');
% Gamma region limits
cbX = linspace(-xLim, xLim, 1e2);
plot(cbX, GammaLim*ones(size(cbX)),'--','color',colors(3,:),'LineWidth',2)
plot(cbX, -GammaLim*ones(size(cbX)),'--','color',colors(3,:),'LineWidth',2)
% Tracking ellipse
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
patchell(P/rhoStar, [iR(idx),vR(idx)], 'patch', colors(4,:));

hold off
axis([-xLim xLim -yLim yLim])
set(gca, 'FontSize',14)

%% Output Voltage
subplot(3,2,2)
plot(t(1:idx), vR(1:idx), 'color', colors(2,:), 'LineWidth', scale*1)
hold on
scatter(t(idx), vR(idx), scale*80, colors(2,:), 'filled');
plot(t(1:idx), vC(1:idx), 'color', colors(1,:), 'LineWidth', scale*1)
scatter(t(idx), vC(idx), scale*50, colors(1,:), 'filled');
hold off
axis([min(t) max(t) min(min(vC),min(vR)) - 10 max(max(vC), max(vR)) + 10])
set(gca, 'FontSize',14)

%% Output Current
subplot(3,2,4)
plot(t(1:idx), iR(1:idx), 'color', colors(2,:), 'LineWidth', scale*1)
hold on
scatter(t(idx), iR(idx), scale*80, colors(2,:), 'filled');
plot(t(1:idx), iL(1:idx), 'color', colors(1,:), 'LineWidth', scale*1)
scatter(t(idx), iL(idx), scale*50, colors(1,:), 'filled');
hold off
axis([min(t) max(t) min(min(iL),min(iR)) - 10 max(max(iL), max(iR)) + 10])
set(gca, 'FontSize',14)

%% Switch variable
subplot(3,2,6)
plot(t(1:idx), q(1:idx), 'color', colors(5,:), 'LineWidth', scale*1)
hold on
scatter(t(idx), q(idx), scale*50, colors(5,:), 'filled');
hold off
axis([min(t) max(t) min(q)-0.1 max(q)+0.1])
set(gca, 'FontSize',14)
%xlabel('$t$ [s]','interpreter','latex')
%ylabel('$q$','interpreter','latex')