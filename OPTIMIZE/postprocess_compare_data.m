function [t_obs,theta_observed]=postprocess_compare_data(inputParams)


% Ks_val=convexCombo(s_param,t_param,3220,108.5,405.1);
% Ks_val=Ks_val*1e-3/(24*3600);
% n_val=convexCombo(s_param,t_param,2.503,1.149,1.649);
% alpha_val=convexCombo(s_param,t_param,3.321,2.711,0.8294);
% eta_val=convexCombo(s_param,t_param,-0.8653,-5.153,0.5452);
% theta_r_val=convexCombo(s_param,t_param,0.0515,0.0961,0.0506);
% theta_s_val=convexCombo(s_param,t_param,0.3769,0.4616,0.5204);

Ks_val=inputParams(1);
n_val=1.4;
alpha_val=inputParams(2);
eta_val=0.5;
theta_r_val=0.08;
theta_s_val=0.4;


params.Ks = Ks_val;
params.n = n_val;
params.alpha=alpha_val;
params.theta_r = theta_r_val;
params.theta_s = theta_s_val;
params.eta     = eta_val;
params.capacity_correction=9.81e-7; % From Mathias et al., under Equation (17)


% load data from file
temp2=load('data.mat');
% temp2=load('johnstown_L1p8m_BC_saturated.mat');
t=temp2.temp.t; % time in days
z=temp2.temp.z; % z coordinate
h_tz=temp2.temp.h_tz; % spacetime array of pressure head, in SI units.
theta_tz=temp2.temp.theta_tz; % spacetime array of theta values, in units of volume/volume.

% *************************************************************************
myTable=readtable('obs_15cm_johnstown.xlsx');

% time = T{:,1};         % datetime column
Tobs = myTable{:,2:end}; % numeric columns

t_obs=Tobs(:,1);
head_hPa=Tobs(:,2);

t_obs=(t_obs/365)+1998; % observations start on Jan 1st 1998

rho_g=1000*9.8;
head_SI=-head_hPa*100/(rho_g);

% *************************************************************************

L=max(z); % length of column.
% Caution!  z=L is the top.
z15=L-0.15;
[~,ix]=min(abs(z-z15));

t1=(t/365)+1997; % met data starts on Jan 1st 1997
[~,ix1]=min(abs(t1-1998));
[~,ix2]=min(abs(t1-2001));
t_lwr=t1(ix1+1);
t_upr=t1(ix2-1);
dt=(t_upr-t_lwr)/365;
t_interp=t_lwr:dt:t_upr;

h_model=interp1(t1(ix1:ix2),h_tz(ix1:ix2,ix),t_interp,'spline');
theta_model=interp1(t1(ix1:ix2),theta_tz(ix1:ix2,ix),t_interp,'spline');

% *************************************************************************

[t_obs_u, idx] = unique(t_obs);
head_SI_u=head_SI(idx);
h_observed=interp1(t_obs_u,head_SI_u,t_interp);

theta_observed=0*t_interp;

for i=1:length(h_observed)
    h_val=h_observed(i);
    [theta_out, ~,~] = vg_model(h_val, params);
    theta_observed(i)=theta_out;
end

% *************************************************************************
% % % % Plotting - theta
% % 
% Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% plot(t_interp,theta_observed,'DisplayName','Observations','Marker','o');
% hold on
% 
% plot(t_interp,theta_model,'DisplayName','Model');
% 
% set(gca,'xtick',[1998 1999 2000 2001])
% 
% grid on
% 
% ylabel('\theta(z=15cm) [m]')
% % xlim([1998 2001])
% 
% box(axes1,'on');
% grid(axes1,'on');
% hold(axes1,'off');
% % Create legend
% legend(axes1,'show');

% *************************************************************************
% % Plotting - h
% 
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(t_interp,h_observed,'DisplayName','Observations','Marker','o');
hold on

plot(t_interp,h_model,'DisplayName','Model');

set(gca,'xtick',[1998 1999 2000 2001])

grid on

ylabel('h(z=15cm) [m]')
% xlim([1998 2001])

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend(axes1,'show');


end