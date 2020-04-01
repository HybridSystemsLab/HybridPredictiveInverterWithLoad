function inside = D_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: D_inverter.m
%
% Description: Jump set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L R Cap Vdc omega rhoStar deltaBar epsTilt lambda varepsilon
%% Control states
q   = x(1);
iL  = x(2);
vC  = x(3);
iR  = x(4);
vR  = x(5);

%% Estimator's states
QQ = x(12);
ell = x(17); % Connectivity of the load
theta = x(18); % Value of the load
counter = x(16); %counter

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


%% Switching rules
inside = 0;
if rhoStar <= V_B && V_B <= deltaBar
    inside = dV_B >= -lambda*k*V_B;
end

if ( (QQ >= varepsilon) && (ell == 1) && (counter <= 1) )
    inside  = 1;
end
 
if ( (V_B <= rhoStar) && (ell == 0) )
   inside  = 1;
end
  
end


