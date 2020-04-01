function VDC = VDCnoise(t)

VDC0 = 220;
% Step changes every 0.1s
VDC = VDC0;
if t >= 0.1 && t < 0.2
    VDC = VDC-20;
elseif t >= 0.2 && t < 0.3
    VDC = VDC+20;
elseif t >= 0.3 && t < 0.4
    VDC = VDC+40;
else
    VDC = VDC;
end

% Sinusoidal noise with bounded amplitude 
A = .02*VDC;
f = 200;
VDC = VDC + A*sin(2*pi*f*t);

end


