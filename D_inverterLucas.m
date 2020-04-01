function inside = D_inverterLucas(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: D_inverter.m
%
% Description: Jump set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L R Cap Vdc omega rhoStar deltaBar epsTilt lambda 
%% states
q   = x(1);
iL  = x(2);
vC  = x(3);
iR  = x(4);
vR  = x(5);
ell = x(7);
theta = x(8);
%% System imtermediate state
ei = iL-iR;
ev = vC-vR;
e = [ei; ev];
P = [1, epsTilt/2*(1-ell); epsTilt/2*(1-ell), (Cap*omega)^2];
A = [-R/L, -omega^2*Cap; 1/Cap, -theta/Cap*ell];
nu = Vdc/L*q  - R/L*iR + (L*Cap*omega^2-1)/L*vC + (-theta*iR + vR*theta^2)/Cap*ell ;
V_B = e'*P*e;
dV_B = e'*(A'*P + P*A)*e + 2*e'*P*[nu;0];
 k = (1-ell) * R/L + 2* ell * min ( R/L , theta/Cap ); % \lambda in the paper


%% Switching rule
if rhoStar <= V_B && V_B <= deltaBar
    inside = dV_B >= -lambda*k*V_B;
else
    inside = 0;
end
