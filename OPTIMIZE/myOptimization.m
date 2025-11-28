function [p_opt,fval]=myOptimization()

% main: Optimize three parameters p1, p2, p3 with bounds.

% parpool   % only needed once per MATLAB session

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

% parameters:
% Silty: s=0, t=0;
% Sandy: s=1, t=0;
% Clay: s=0, t=1;

lb=[1e-6,0.7];
ub=[5e-6,1];

p0=(lb+ub)/2;

% Clear previous log
logFile = 'cost_log.txt';
if exist(logFile,'file')
    delete(logFile);
end

% *************************************************************************
% Setting up intermediate parameters by calling richards_pdepe evaluated at
% paramer values p0.

disp('Setting up intermediate parameters')
% [z,t,h_tz,theta_tz,sol]=richards_pdepe(p0);

% load data from file
temp2=load('data.mat');
% temp2=load('johnstown_L1p8m_BC_saturated.mat');
t=temp2.temp.t; % time in days
z=temp2.temp.z; % z coordinate
h_tz=temp2.temp.h_tz; % spacetime array of pressure head, in SI units.
theta_tz=temp2.temp.theta_tz; % spacetime array of theta values, in units of volume/volume.


L=max(z); % length of column.
% Caution!  z=L is the top.
z15=L-0.15;
[~,ix]=min(abs(z-z15));

t1=(t/365)+1997; % met data starts on Jan 1st 1997
[~,ix1]=min(abs(t1-1998));
[~,ix2]=min(abs(t1-2000));
t_lwr=t1(ix1+1);
t_upr=t1(ix2-1);
dt=(t_upr-t_lwr)/365;
t_interp=t_lwr:dt:t_upr;

h_model=interp1(t1(ix1:ix2),h_tz(ix1:ix2,ix),t_interp,'spline');

[t_obs_u, idx] = unique(t_obs);
head_SI_u=head_SI(idx);
h_observed=interp1(t_obs_u,head_SI_u,t_interp);

disp('done')

J_init = sum( (h_observed-h_model).^2)/length(t_interp);

display(strcat('preliminary value of cost function is J=',num2str(J_init)))

% *************************************************************************


% % % Optimization options (optional)
% % opts = optimoptions('fmincon', ...
% %     'Display', 'iter', ...
% %     'Algorithm', 'sqp');
% % 
% % % Run optimization
% % [p_opt, fval] = fmincon(@my_cost, p0, [], [], [], [], lb, ub, [], opts);

% Options for surrogate optimization
opts = optimoptions('surrogateopt', ...
    'UseParallel', true, ...
    'Display', 'iter', ...         % Show iteration info
    'MaxFunctionEvaluations', 100); % Limit expensive PDE calls for demo

% Run surrogate optimization
[p_opt, fval] = surrogateopt(@my_cost, lb, ub, opts);

fprintf('Optimal parameters:\n');
disp(p_opt)
fprintf('Optimal cost: %.6f\n', fval);

% *************************************************************************
% my_cost: User-defined cost function
% p = [p1; p2; p3]

    function J = my_cost(p)

        % p1 = p(1);
        % p2 = p(2);
        % p3 = p(3);

        fid = fopen('cost_log.txt','a');
        % fprintf(fid, 'Evaluating p1, p2 = %f %f\n', p(1),p(2));
        % fprintf(fid, 'Evaluating p= %f\n', p);


        [z_local,t_local,h_tz_local,theta_tz_local,sol_local]=richards_pdepe(p);
        
        h_model_local=interp1(t1(ix1:ix2),h_tz_local(ix1:ix2,ix),t_interp,'spline');
        
        J = sum( (h_observed-h_model_local).^2)/length(t_interp);
        fprintf(fid, 'cost J= %f\n', J);
        fclose(fid);
    end

% *************************************************************************

end
