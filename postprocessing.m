%% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: postprocessing.m
%
% Description: processing and plot of simulation data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Cap omega ePerc H Vdc GammaLim rhoStar deltaBar epsTilt
colors = lines;

q = x(:,1);
iL = x(:,2);
vC = x(:,3);
iR = x(:,4);
vR = x(:,5);
tau = x(:,6);

%% Plot states separately
figure
subplot(311)
plot(t, q, 'color', colors(1,:));
xlabel('Time (s)')
ylabel('q')
title('Switching variable')
axis([0 t(end) -Inf Inf])

subplot(312)
plot(t, iL, 'color', colors(2,:));
xlabel('Time (s)')
ylabel('iL [A]')
title('Inductor Current')
axis([0 t(end) -Inf Inf])

subplot(313)
plot(t, vC,t, GammaLim* ones (size(t)), 'color', colors(3,:));
xlabel('Time (s)')
ylabel('vC [V]')
title('Output Voltage')
axis([0 t(end) -Inf Inf])

box on
print -depsc2 states.eps

%% Simulation Video
VideoScript;

%% Phase Plane - Final
figure('Units', 'inches', 'Position', [10 10 10 8])
hold on
xLim = 1.1*GammaLim*Cap*omega;
yLim = 1.1*GammaLim;

pidx = 1:10:length(t);
% Tracking ellipse at t = end
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
h(1) = patchell(P/deltaBar, [iR(idx),vR(idx)], 'patch', colors(5,:));
% Inverter state
h(2) = plot(iL(pidx), vC(pidx),'color',colors(1,:),'LineWidth', 1.5);
% Reference state
h(3) = plot(iR(pidx), vR(pidx),'color',colors(2,:),'LineWidth', 1.5);
% Boundaries of the Gamma region
cbX = linspace(-100, 100, 1e2);
plot(cbX, GammaLim*ones(size(cbX)), '--', 'color', colors(4,:))
plot(cbX, -GammaLim*ones(size(cbX)), '--', 'color', colors(2,:))

legend(h(1:3), {'\qquad\qquad', '',''},'interpreter','latex','FontSize',14)
box on
set(gca, 'FontSize',14)
axis([-xLim xLim -yLim yLim])
print -depsc2 LocalControllerPP.eps

%% Error Phase Plane
figure('Units', 'inches', 'Position', [10 10 10 8])
hold on
colors = lines;
pidx = 1:length(t);
% Tracking ellipse
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
patchell(P/rhoStar, [0,0], '', colors(4,:));
% Estimate of the region of attraction
P = [H, epsTilt/2; epsTilt/2, (Cap*omega)^2];
patchell(P/deltaBar, [0,0], '', colors(5,:));
% Error state
h(1) = plot(iL(pidx)-iR(pidx), vC(pidx)-vR(pidx),...
    'color',colors(1,:),'LineWidth', 1.5);
% Error initial condition
scatter(iL(1)-iR(1), vC(1)-vR(1), 150,colors(1,:), 'filled')
set(gca, 'FontSize',14)
box on, grid on
%print -depsc2 LocalControllerError.eps

%% Spectral analysis
% Fast Fourier Transform
[f, vcFFT] = nfftAlgorithm(t, vC, solver);
figure;
plot(f, vcFFT);
title('Fourier transform');
xlabel('frequency [Hz]');
ylabel('vC(f) [dB]');
grid on;
axis([0 500 -Inf Inf])

% Total Harmonic Distortion
figure
thd(vC, length(vC)/t(end));
[thd_db,~,~] = thd(vC, length(vC)/t(end));
percent_thd = 100*(10^(thd_db/20));
disp(percent_thd);
axis([0 1 -Inf Inf])