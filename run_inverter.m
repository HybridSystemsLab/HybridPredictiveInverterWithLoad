%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file
% Project: Hybrid Predictive Inverter
%
% Filename: run_inverter
%
% Description: simulation of a Single-Phase DC/AC inverter controlled by 
%              a Hybrid Predictive Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
figure 
for zz = 1:1:2
    for zzz = 1:1:2

%% Useful Constants
global omega R L Cap Vdc a b T rhoStar NumErr RL H Vout RL_real
global VoutMax deltaBar Levels lgc epsTilt GammaLim lambda varepsilon
global vc1jump il1jump vr1jump ir1jump vc2jump il2jump vr2jump ir2jump
%% Circuit parameters
varepsilon = 0.4;
R = 1;                  % resistance ||Ohms||
L = 2e-3;               % inductance ||Henry's||
Cap = 1.063e-3;         % capacitance ||Farads||
Vdc = 220;              % DC input voltage ||Volt||
f = 50;                 % frequency ||Hz||
omega = f*2*pi;         % angular freqency
    RL = 60;           % best guess 
    RL_real = 100;      % load resistance ||Ohms|| the real one

    
    Levels = 3;             % inverter levels (see multilevel inverters)
% Tracking ellipse parameters
Vout = 100;             % desired output voltage ||Volt||
b = Vout;               % reference ellipse: vB axis
a = b*Cap*omega;        % reference ellipse: iB axis
ePerc = .1;             % tracking error percentage
rhoStar = 4; % delta in the paper 
H = 1;                  % aspect ratio of the tracking ellipse
epsTilt = R*Cap/L;      % mixed term of function V
T = 0.001;             % prediction horizon ||s||
lgc = 1;                % used for the impulsive resets of z_r (see G.m)

%% Controller constraints (to be modified)
lambda = .1;
k = abs(L*Cap*omega^2-1);
VoutMax = (Vdc/k-sqrt(rhoStar/((Cap*omega)^2 ...
    - (.5*epsTilt)^2)))*(k/(k + omega*R*Cap));
deltaBar = ((Cap*omega)^2-(.5*epsTilt)^2)*...
    ((Vdc-Vout*(omega*R*Cap + k))/k)^2;
GammaLim = (Vdc - omega*R*Cap*Vout)/k;
if Vout > VoutMax
    error('Vout <= Vout_Max: constraint not satisfied')
end
if deltaBar < rhoStar
    error('rhoStar <= deltaBar: constraint not satisfied')
end
if lambda <= 0 || lambda >= 1
    error('lambda \in (0,1): constraint not satisfied')
end
NumErr = 1e-10;

%% Initial setting
global plotSolutions
plotSolutions = 0;  % Plot predicted trajectories at each jump
                    % (requires setting a breakpoint in g_inverter.m
if plotSolutions
    figure('Units', 'inches', 'Position',[10 10 12 10]), hold on;
end
global PredictionMethod
% FixedHorizon:     predicts solutions for all t \in [0, T]
% EventDetection:   predicts solutions until they hit the jump set
PredictionMethod = 'EventDetection';
global simulateVDCnoise simulateRefResets
simulateVDCnoise = 0;   % Adds two noise signals to the input voltage VDC 
simulateRefResets = 0;  % Simulates impulsive resets of the reference zr

aEll = sqrt(deltaBar)/sqrt(H);
bEll = sqrt(deltaBar)/(Cap*omega);
thetaEll = 1/2*atan(epsTilt/(H-(Cap*omega)^2));
th = 0;
thR = 0;
rho0 = rhoStar;

%% Initial conditions (control, circuit, reference)
ell0= 0;
q0 = 0;
iR0 = omega*Cap*Vout*cos(thR)  + Vout * sin(thR) * ell0 / RL ;  
vR0 = Vout*sin(thR);
iL0 = iR0 + (-1)^zz * 4; % initial current 
vC0 = vR0 + (-1)^zzz * 5; % initial voltage
%iL0 = iR0 ; % initial current 
%vC0 = vR0 - 6; % initial voltage
tau0 = 0;
x0 = [q0; iL0; vC0; iR0; vR0; tau0];
%% Initial conditions (estimator)
xo0     = [iL0 vC0]';   % observer's state hatz0 in the paper 
th0     = 1/RL;         % parameter's estimate hatthetae0 in the paper 
LL0     = [0.1 0.1]';   % auxiliary state L
QQ0     = 0;            % auxiliary state Q
eta0    = [0.1 0.2]';   % auxiliary state \eta
gamma0  = 0.1;          % auxiliary state \Gamma 
y0 = [xo0;th0;LL0;QQ0;eta0;gamma0]; % x_e in the paper

%% Initial conditions (connect load, counter, load estimate)
counter0 = 0; % counter = {0, 1, 2}, r in the paper
theta0 = 1/RL; % hattheta0 according to the paper

%% Simulation
TSPAN = [0 .06];
JSPAN = [0 1e4];
rule = 1;
solver = 'ode45';
options = odeset('RelTol',1e-5,'MaxStep',1e-5);
[t, j, x] = HyEQsolver(@f_inverter, @g_inverter, @C_inverter, @D_inverter, ...
   [x0;y0;counter0; ell0; theta0],TSPAN,JSPAN,rule,options, solver);


%% The plot of e_V function of e_i:
iL = x(:,2);
vC = x(:,3);
iR = x(:,4);
vR = x(:,5);
ei_post_hat = iL - iR;
ev_Post_hat = vC - vR;

%figure 
hp = plot(ei_post_hat, ev_Post_hat,'LineWidth',1);
grid on;
hold on
xlabel('$e_i [A]$','interpreter','latex', 'fontsize',20);
ylabel('$e_v [V]$','interpreter','latex', 'fontsize',20)

f=@(x,y) x.^2 + (Cap * omega).^2 * y.^2 + epsTilt .* x .* y; 
[X,Y]=meshgrid(-10:0.01:10);
z=f(X,Y);
contour(X,Y,z,[rhoStar,rhoStar],'k')
hold on
f=@(x,y) x.^2 + (Cap * omega).^2 * y.^2;
[X,Y]=meshgrid(-10:0.01:10);
z=f(X,Y);
contour(X,Y,z,[rhoStar,rhoStar],'r')
hold on

ei1jump = il1jump - ir1jump;
ev1jump = vc1jump - vr1jump;   
ei2jump = il2jump - ir2jump;
ev2jump = vc2jump - vr2jump; 

if ( zz == 2 && zzz == 2)
    plot(ei1jump,ev1jump,'Marker','o','MarkerSize',10,'MarkerFaceColor','m');
    hold on 
    plot(ei2jump,ev2jump,'Marker','o','MarkerSize',10,'MarkerFaceColor','c');
    hold on
end

axis([-8 8 -10 10]);

%% Plot solution
% postprocessing_Haoyue;
   end
end
 

