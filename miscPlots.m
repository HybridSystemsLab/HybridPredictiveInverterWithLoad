%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: miscPlots.m
%
% Description: plots of two robustness tests for the closed-loop system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VDC noise
% Profile plot of the perturbed VDC signal and the resulting vC and vR.
figure('Units', 'inches', 'Position', [10 10 10 3])
hold on
plotidx = 1:40:length(t);
plot(t(plotidx), vC(plotidx),'color',colors(1,:),'LineWidth', 1.5);
plot(t(plotidx), vB(plotidx),'color',colors(2,:),'LineWidth', 1);
legend( {'\qquad\qquad', '',''}, 'interpreter','latex', 'FontSize', 14)
set(gca, 'FontSize',14)
ylim([-120,120])
box on
print -depsc -r300 VdcnoiseStates.eps

figure('Units', 'inches', 'Position', [10 10 10 3])
plot(t(plotidx), arrayfun(@VDCnoise, t(plotidx)), 'color', colors(5,:),'LineWidth', 1.5)
set(gca, 'FontSize',14)
gcf.PaperPosition = [0 0 10 5];
box on
print('Vdcnoise.eps', '-depsc2', '-r300')

%% Reference resets
% Plot of vC and vR subjected to impulsive resets.
figure('Units', 'inches', 'Position', [10 10 10 5])
hold on
plotidx = 1:40:length(t);
hold on
h(5) = plot(t(plotidx), vC(plotidx),'color',colors(1,:),'LineWidth', 1.5);
h(6) = plot(t(plotidx), vB(plotidx),'color',colors(2,:),'LineWidth', 1);
legend(h(5:6), {'\qquad\qquad', ''}, 'interpreter','latex', 'FontSize', 14)
set(gca, 'FontSize',14)
print -depsc -r300 referenceResets.eps
grid on
box on