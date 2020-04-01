function xdot = f_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: f_inverter.m
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Vdc L R Cap omega simulateVDCnoise RL_real Vout
%% control states
q   = x(1);
iL  = x(2);
vC  = x(3);
iR  = x(4);
vR  = x(5);
tau = x(6);

%% Estimator's states
xo    = x(7:8); % iL_hat  vC_hat
th    = x(9);  % hatthetae
LL    = x(10:11);  
eta   = x(13:14);  
ell = x(17); % Connectivity of the load
theta = x(18); % Value of the load hattheta


%% VDC noise 
if simulateVDCnoise
    Vdc = VDCnoise(tau);  
end

%% control dynamics
tt = 1/RL_real;
qdot = 0;
iLdot = Vdc/L*q -R/L*iL - 1/L*vC;
vCdot = (1/Cap)*iL - (vC*tt/Cap)*ell;
iRdot = -omega^2*Cap*vR + (theta*iR-vR*theta^2)*ell/Cap ;
vRdot = (1/Cap)*iR - (vR/Cap)*theta*ell;
taudot = 1;
xcdot = [qdot; iLdot; vCdot; iRdot; vRdot; taudot];

%% estimator dynamics
AA = [-R/L -1/L; 1/Cap 0];
xp    = [iL;vC];
gx    = [0; -1/Cap*vC];
KK=1;
m = [Vout/L*q; 0];
fx = AA*xp;
h = 0;
xodot = fx + m+  gx*th + KK*(xp - xo) + LL*h;
thdot = h;
LLdot = gx - KK*LL;
QQdot = LL'*LL;
etadot  =  - KK*eta;
gammadot=  LL'*(LL*th + xp - xo - eta);
xedot = [xodot; thdot; LLdot; QQdot; etadot; gammadot];      

%% extra states dynamics
counterdot = 0;
elldot = 0;
thetadot = 0;


%% overall dynamics
xdot = [xcdot; xedot; counterdot; elldot; thetadot];

end