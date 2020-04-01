%% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Hybrid Predictive Inverter
%
% Name: postprocessing.m
%
% Description: processing and plot of simulation data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rhoStar  Cap omega epsTilt L Vdc R Vout
% close all; 
colors = lines;
q = x(:,1);
iL = x(:,2);
vC = x(:,3);
iR = x(:,4);
vR = x(:,5);
ei_post_hat = iL - iR;
ev_Post_hat = vC - vR; 
ell = x(:,17);
th = x(:,9);
theta = x(:,18);
counter = x(:,16);

%% The plot of e_V function of e_i:
% figure 
% hp = plot(ei_post_hat, ev_Post_hat,'LineWidth',1);
% grid on;
% hold on
% xlabel('$e_i [A]$','interpreter','latex', 'fontsize',20);
% ylabel('$e_v [V]$','interpreter','latex', 'fontsize',20)
% axis([-8 8 -15 15]);
% grid on
% f=@(x,y) x.^2 + (Cap * omega).^2 * y.^2 + epsTilt .* x .* y; 
% [X,Y]=meshgrid(-8:0.01:8);
% z=f(X,Y);
% contour(X,Y,z,[rhoStar,rhoStar])
% hold on
% f=@(x,y) x.^2 + (Cap * omega).^2 * y.^2;
% [X,Y]=meshgrid(-8:0.01:8);
% z=f(X,Y);
% contour(X,Y,z,[rhoStar,rhoStar],'g')
% hold on

%% V versus delta
e_t = [ei_post_hat ev_Post_hat];
v_post = zeros(1,length(e_t));
rho = rhoStar * ones(1,length(t));
for i=1:length(e_t)
    P = [1, epsTilt/2*(1-ell(i)); epsTilt/2*(1-ell(i)), (Cap*omega)^2];
    v_post(i) =  e_t(i,:)*P*e_t(i,:)'; 
end

figure 
hp = plot(t, v_post,t, rho,'LineWidth',1);
legend(hp((2)), '$\delta$','interpreter','latex', 'fontsize',20);
hold off;
grid on;
xlabel('t[sec]','interpreter','latex', 'fontsize',20);
ylabel('$V(e(t))$','interpreter','latex', 'fontsize',20);

%% Estimator's Plot
figure 
plotflows(t,j,th);
grid on;
xlabel('t[sec]','interpreter','latex', 'fontsize',20);
ylabel('$\hat{\theta}_{e}$ [$\Omega$]','interpreter','latex', 'fontsize',20);
axis([0 0.004 0 0.3] );

%% The plot of thetahate and theta
% figure 
% hp = plot(t, th,t,theta,'LineWidth',1);
% grid on;
% xlabel('t[sec]','interpreter','latex', 'fontsize',20);
% ylabel('$\hat{\theta}_e(t)$','interpreter','latex', 'fontsize',20)

%% The plot of ell
% figure 
% hp = plot(t, ell,'LineWidth',1);
% grid on;
% xlabel('t[sec]','interpreter','latex', 'fontsize',20);
% ylabel('$\ell(t)$','interpreter','latex', 'fontsize',20)

%% Plot showing if vc is in Gamma
% k = abs(L*Cap*omega^2-1);
% Gamma = ( (Vdc - R * Cap * omega * Vout)/k )* ones( size(vC) );
% figure 
% hp = plot(t, abs(vC),t, Gamma,'LineWidth',1);
% grid on;
% xlabel('t[sec]','interpreter','latex', 'fontsize',20);
% ylabel('$v_c(t)$','interpreter','latex', 'fontsize',20);

%% plot q function of time 
%figure 
%hp = plot(t, q,'LineWidth',1);
%grid on;
%xlabel('t[sec]','interpreter','latex', 'fontsize',20);
%ylabel('$q(t)$','interpreter','latex', 'fontsize',20);

%% Total Harmonic Distortion  for VC %figure 6
% figure
% thd(vC, length(vC)/t(end));
[thd_db,~,~] = thd(vC, length(vC)/t(end));
percent_thd = 100*(10^(thd_db/20));
disp('The V_C thd% is: ')
disp(percent_thd);
% axis([0 1 -Inf Inf])

%% Total Harmonic Distortion for IL %figure 6
% figure
% thd(iL, length(iL)/t(end));
[thd_db_iL,~,~] = thd(iL, length(iL)/t(end));
percent_thd_iL = 100*(10^(thd_db_iL/20));
disp('The i_L thd% is: ')
disp(percent_thd_iL);
% axis([0 1 -Inf Inf])

%% 
disp('The amount of switches: ');
disp(size( find ( diff (q) ~=0 ) ));