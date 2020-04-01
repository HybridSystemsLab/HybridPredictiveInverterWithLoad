%% Simulation of the closed-loop system for frequency analys

lvls = [3,5,7,9,11];        % Array with number of levels
nLev = length(lvls);
nSim = 8;                   % Number of simulations
testCell = cell(l,nSim);    % Cell array containing simulation data
TSPAN = [0 5];
for l = 1:nLev
    for i = 1:nSim
        iR0 = omega*Cap*Vout*cos((i-1)/2*pi);
        vR0 = Vout*sin((i-1)/2*pi);
        iL0 = iR0;
        vC0 = vR0;
        x0 = [q0; iL0; vC0; iR0; vR0; tau0];
        Levels = lvls(l);
        [t, j, x] = HyEQsolver(@f_inverter, @g_inverter, @C_inverter,....
            @D_inverter, x0, TSPAN, JSPAN, rule, options, solver);
        % Fill the cell array
        testCell{l,i}.t = t;
        testCell{l,i}.q = x(:,1);
        testCell{l,i}.iL = x(:,2);
        testCell{l,i}.vC = x(:,3);
    end
end

%% FFT computation and plots
figure('Units', 'inches', 'Position', [10 10 10 3])
hold on
colors = lines;
fFFT = cell(1,3);
sFFT = cell(1,3);
for l = 1:nLev
    vFFT = 0;
    f = 0;
    for i = 1:nSim
        % Fast Fourier Transform
        time = testCell{l,i}.t;
        signal = testCell{l,i}.vC;
        [f, vcFFT] = nfftAlgorithm(time, signal, solver);
        plot(f, 20*log10(vcFFT), 'LineWidth',1.5,'color',colors(i,:));
        box on
        axis([0 350 -70 50])
    end
end

set(gca, 'FontSize',14)
print -depsc2 FFT.eps

%% Total Harmonic Distortion
% Initialize arrays for data on THD
percentTHDvoltage = zeros(nLev, nSim);
percentTHDcurrent = zeros(nLev, nSim);
for l = 1:nLev
    for i = 1:nSim
        % Fast Fourier Transform
        time = testCell{l,i}.t;
        signal = testCell{l,i}.vC;
        [thd_db,~,~] = thd(signal, length(signal)/time(end));
        percentTHDvoltage(l,i) = 100*(10^(thd_db/20));
        signal = testCell{l,i}.iL;
        [thd_db,harmpow,harmfreq] = thd(signal, length(signal)/time(end));
        percentTHDcurrent(l,i) = 100*(10^(thd_db/20));
    end
end

disp(percentTHDvoltage);
disp(percentTHDcurrent);