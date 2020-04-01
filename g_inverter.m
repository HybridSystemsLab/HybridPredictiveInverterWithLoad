function xplus = g_inverter(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Name: g_inverter.m
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plotSolutions PredictionMethod simulateRefResets epsTilt rhoStar
global T omega Cap R L Vdc Levels lgc Vout  varepsilon lambda
global vc1jump il1jump vr1jump ir1jump vc2jump il2jump vr2jump ir2jump
%% states
q = x(1);
iL = x(2);
vC = x(3);
iR = x(4);
vR = x(5);
tau = x(6);

% Estimator's states
xo    = x(7:8);
th    = x(9); % hatthetae
LL    = x(10:11);
QQ    = x(12);
eta   = x(13:14);
gamma = x(15);
counter = x(16);  %counter
ell = x(17); % Connectivity of the load
theta = x(18); % hattheta

%% safety jump
%% initilize jump value
xp = [iL;vC];
iLplus = iL; 
qplus = q;
vCplus = vC; 
xoplus = xo; 
thplus = th; 
LLplus = LL; 
QQplus = QQ; 
etaplus= eta; 
gammaplus = gamma; 
counterplus = counter;
ellplus = ell; 
thetaplus = theta; 


%% Jump map for the controller qplus
ei = iL-iR;
ev = vC-vR;
e = [ei; ev];
P = [1, epsTilt/2*(1-ell); epsTilt/2*(1-ell), (Cap*omega)^2];
A = [-R/L, -omega^2*Cap; 1/Cap, -theta/Cap*ell];
nu = Vdc/L*q  - R/L*iR + (L*Cap*omega^2-1)/L*vC + (-theta*iR + vR*theta^2)/Cap*ell ;
V_B = e'*P*e;
dV_B = e'*(A'*P + P*A)*e + 2*e'*P*[nu;0];
 k = (1-ell) * R/L + 2* ell * min ( R/L , theta/Cap ); % \lambda in the paper

 if (rhoStar <= V_B) && (dV_B >= -lambda*k*V_B)
TPspan = [0 T];
q0 = linspace(-1,1,Levels);
Tif = zeros(size(q0)); % Initialize time to impact function
% Compute qhat according to the map Ghat
qbar = ( R*iR - (L*Cap*omega^2-1)*vC - L * (vR * theta^2 - iR * theta ) * ell / Cap ) /Vdc;
if ei < -(R*Cap/(2*L)*ev)*(1-ell)
    qhat = q0(q0 >= qbar);
elseif ei > -(R*Cap/(2*L)*ev)*(1-ell)
    qhat = q0(q0 <= qbar);
else
    qhat = q0;
end

switch PredictionMethod
    case 'FixedHorizon'
        options = odeset('RelTol',1e-5,'MaxStep',1e-4,'Events', []);
        if plotSolutions
            qhat = q0;
        end
        for i = 1:length(qhat)
            x0 = [qhat(i), iL, vC, iR, vR, ell, theta];
            [tP, xP] = ode45(@(t, xP) f_inverterLucas(xP), TPspan, x0, options);
            Tif(i) = 1;
            % Increment Iq if xP \in C
            while Tif(i) <= size(xP,1) && ~D_inverterLucas(xP(Tif(i),:))
                Tif(i) = Tif(i) + 1;
            end
            if plotSolutions
                len(i) = Tif(i);
                solPlot{i}.t = tP;
                solPlot{i}.x = xP;
            end
        end
        qplus = qhat(Tif == max(Tif));
        
    case 'EventDetection'
        options = odeset('RelTol',1e-5,'MaxStep',1e-4,'Events', ...
            @odeJumpEvent);
        if plotSolutions
            qhat = q0;
        end
        for i = 1:length(qhat)
            x0 = [qhat(i), iL, vC, iR, vR, 0, ell, theta];
            if ~D_inverterLucas(x0) 
                [tP, xP]=ode45(@(t, xP) f_inverterLucas(xP),TPspan,x0,options);
                Tif(i) = tP(end);
            else
                tP = 0;
                xP = x0;
            end
            if plotSolutions
                len(i) = length(tP);
                solPlot{i}.t = tP;
                solPlot{i}.x = xP;
            end
        end
        qplus = qhat(Tif == max(Tif));
    otherwise
        qplus = q;
end

% If the jump map is set-valued, randomly select one qplus among G(x)
if length(qplus) > 1
    qplus = qplus(ceil(rand*length(qplus)));
end


%% Plot predicted trajectories
if plotSolutions
    plotPredTrajectories;
    % (set breakpoint here)
end

 end
 
 
 
 
 
%% Simulate impulsive reference resets
if simulateRefResets
    if tau >= 0.05
        Vout = Vout - lgc*Vout/2;
        lgc = -lgc;
        theta = rand*2*pi;
        iRplus = Cap*omega*Vout*cos(theta) + Vout * sin(theta) * ell * theta ;
        vRplus = Vout*sin(theta);
        tauplus = 0;
    else
        tauplus = tau;
        iRplus = iR;
        vRplus = vR;
    end
else
    iLplus = iL;
    vCplus = vC;
    iRplus = iR;
    vRplus = vR;
    tauplus = tau;
end


%% Jump for the estimator
 if ( ( QQ >= varepsilon ) && ( ell == 1 ) ) %% Estimator's Jump Map
    xoplus    = xp;
    thplus    = (1/QQ)*gamma;
    LLplus    = [0 0]';
    QQplus    = 0;
    etaplus   = [0 0]';
    gammaplus = 0;
    counterplus = min(counter + 1,2);
    if ( counter == 1 ) %% Record estimator's 2nd jump value for plotting
    vc2jump = vCplus;
    il2jump = iLplus;
    vr2jump = vRplus;
    ir2jump = iRplus;
    thetaplus = thplus;
    end
 end

 %% Jump to switch the load
 
if ( (V_B < rhoStar) && ( ell == 0 ) )
    ellplus = ell + 1;
    vc1jump = vCplus;
    il1jump = iLplus;
    vr1jump = vRplus;
    ir1jump = iRplus;
end

%% finally overall jump
xplus = [qplus; iLplus; vCplus; iRplus; vRplus; tauplus; ...
     xoplus; thplus; LLplus; QQplus; etaplus; gammaplus; ...
    counterplus; ellplus; thetaplus];


end