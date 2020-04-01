function xdot = f_reference(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: f_reference.m
%
% Description: reference dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Cap omega 
% states
iB = x(1);
vB = x(2);

% flow map
iBdot = -Cap*omega^2*vB;
vBdot = 1/Cap*iB;

xdot = [iBdot; vBdot];
end