function [f_vec, FFT_vec] = nfftAlgorithm(t, y, solver)

if ~strcmp(solver, 'ode23')
    % Create a strictly monotonic time vector
    tF = unique(t);
    yF = zeros(size(tF));
    % Reshape the input vector to match the time vector
    k = 0;
    for i = 2:length(t)+1
        if i == length(t)+1
            yF(k) = y(i-1);
        else
            if t(i) ~= t(i-1)
                k = k+1;
                yF(k) = y(i-1);
            end
        end
    end
    % Create the query grid
    Ngrid = length(tF);
    gr = 0:(tF(end)/Ngrid):tF(end);
    % Interpolate the input vector on the query grid
    yInterp = interp1(tF, yF, gr(1:end-1),'spline');
else
    yInterp = y;
    tF = t;
end
% Compute the FFT algorithm on the interpolated signal
FFT=fft(yInterp);
N = length(FFT);
T = tF(end)-tF(1);
f0=1/T;
% Build frequency vector
fr=(0:f0:(N-1)*f0);
% Pick half of the values (FFT symmetry)
FFT_vec = FFT(1:floor(N/2)+1);
f_vec=fr(1:floor(N/2)+1);
% Multiply the amplitude of non-DC components to obtain the same
% spectral density
FFT_vec(2:end) = 2*FFT_vec(2:end);
% Divide by the number of samples
FFT_vec = abs(FFT_vec)/N;

end