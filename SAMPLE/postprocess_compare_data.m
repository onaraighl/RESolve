function [t_obs,theta_obs]=postprocess_compare_data()

% % Data from Johnstown Castle
% 
params.theta_r = 	0.077;
params.theta_s = 0.396;
params.alpha   = 0.894;     % 1/m
params.n       = 1.424;
params.Ks      = 0.195/(24*3600);  % m/s; NOTE - this is the value from Rosetta.
params.eta     = 0.5;
params.capacity_correction=9.81e-7; % From Mathias et al., under Equation (17)

% load data from file
temp2=load('johnstown.mat');
% temp2=load('johnstown_L1p8m_BC_saturated.mat');
t=temp2.temp.t; % time in days
z=temp2.temp.z; % z coordinate
h_tz=temp2.temp.h_tz; % spacetime array of pressure head, in SI units.
theta_tz=temp2.temp.theta_tz; % spacetime array of theta values, in units of volume/volume.

myTable=readtable('obs_15cm_johnstown.xlsx  ');

% time = T{:,1};         % datetime column
Tobs = myTable{:,2:end}; % numeric columns

t_obs=Tobs(:,1);
head_hPa=Tobs(:,2);

rho_g=1000*9.8;
head_SI=-head_hPa*100/(rho_g); % multiply by minus 1 for the correct 
                               % convention whereby h<0 corresponds 
                               % to unsaturated soils.
theta_obs=0*head_SI;
head_obs=head_SI;

for i=1:length(head_SI)
    h_val=head_SI(i);
    [theta_out, ~,~] = vg_model(h_val, params);
    theta_obs(i)=theta_out;
end

L=max(z); % length of column.
% Caution!  z=L is the top.
z15=L-0.15;
[~,ix]=min(abs(z-z15));

% *************************************************************************
% % Plotting - theta
% 
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

t_obs=(t_obs/365)+1998; % observations start on Jan 1st 1998
plot(t_obs,theta_obs,'DisplayName','Observations','Marker','o');
hold on

t=(t/365)+1995; % met data starts on Jan 1st 1995

plot(t,theta_tz(:,ix),'DisplayName','Model');

set(gca,'xtick',[1998 1999 2000 2001])

grid on

ylabel('\theta(depth=15cm)')
xlim([1998 2001])
% ylim([0.382 0.4])

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend(axes1,'show');

% *************************************************************************
% % Plotting - h
% 
% Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% t_obs=(t_obs/365)+1998; % observations start on Jan 1st 1998
% plot(t_obs,head_obs,'DisplayName','Observations','Marker','o');
% hold on
% 
% t=(t/365)+1995; % met data starts on Jan 1st 1995
% 
% plot(t,h_tz(:,ix),'DisplayName','Model');
% 
% set(gca,'xtick',[1998 1999 2000 2001])
% 
% grid on
% 
% ylabel('h(z=15cm) [m]')
% xlim([1998 2001])
% 
% box(axes1,'on');
% grid(axes1,'on');
% hold(axes1,'off');
% % Create legend
% legend(axes1,'show');


end