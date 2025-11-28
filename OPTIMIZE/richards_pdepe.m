function [z,t_vec,h_tz,theta_tz,sol]=richards_pdepe(inputParams)

% function [z,t_vec,h_tz,theta_tz,sol]=richards_pdepe()

% *************************************************************************
% PDE solver to solve the 1D richards equation
%
% (d\theta/dt)=(d/dz)[K(h)( (dh/dz)+1)], z\in [0,L], t>0.     (1)
%
% Lennon Ó Náraigh, with some help from ChatGPT.
% November 2025
% 
% Inputs: null
%
% Outputs:
%  z        - z-coordinate
%  t        - vector of time values
%  h_tz     - spacetime array of values of pressure head (negative for
%             unstaturated).
%  theta_tz - spacetime array of values of theta, in units of mm^3/mm^3.
%  sol      - matlab output (included here as diagnostic).
%
% SI UNITS USED THROUGHOUT
%
% *************************************************************************
% In more detail:
%
% To specifiy the boundary conditions we idenity the flux
% f=K(h)( (dh/dz)+1).
% 
% So the boundary conditions are:
%
% f=[Infiltation flux] at the TOP (top is z=L)
% f=K(h) at the BOTTOM (bottom is z=0).
%
% Initial conditions are provided at t=0, on h.
%
% To solve the Richards Equation, h is identified as the prognostic
% variable.  Equation (1) is transformed to read:
%
% c(h)(dh/dt)=(d/dz)[K(h)( (dh/dz)+1)], z\in [0,L], t>0.      (2)
%
% The capacity c(h) is obtained from the Van Genuchten formula by
% analytical differentiation., c(h)=d\theta(h)/dh.
% 
% Equation (2) is in a standard form which can be handled by Matlab's pdepe
% function:
%
% c(x,t,u)(du/dt)=(d/dx)(f)+s(x,t,u,u_x),
% where f=f(x,t,u,u_x) is the flux function, so we have:
%
% t=t
% x=z
% u=h
% u_x=dh/dz
% f=K(h)( (dh/dz)+1), etc.

% *************************************************************************
% Data from Johnstown Castle, prepared by SK
% Day 1 is Jan 1st 1997, hence myTable{733:end,2:end} below.

myTable=readtable('met_data_johnstown.xlsx');

% time = T{:,1};          % datetime column
Tdata = myTable{733:end,2:end}; % numeric columns, start at Jan 1st 1997

colnumber_EP=3; % column number of evapo-transpiration
colnumber_rainfall=4; % column number of rainfall


% *************************************************************************
% VG parameters:

% Input params
%
% Silty: s=0, t=0;
% Sandy: s=1, t=0;
% Clay: s=0, t=1;

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

params

% params.theta_r = 0.077;
% params.theta_s = 0.396;
% params.alpha   = 0.894;     % 1/m
% params.n       = 1.424;
% params.Ks      = 0.0475/(24*3600);% 0.195/(24*3600);  % m/s; NOTE - this is the value from Rosetta.... =  2.2569e-06 m/s
% params.eta     = 0.5;
% params.capacity_correction=9.81e-7; % From Mathias et al., under Equation (17)

Se_target=0.999;
h_ref = vg_inv_Se(Se_target, params);  % reference level, very close to saturation.

% *************************************************************************
% Plant stress function parameters:

% Note!  Plant root distributin parameters are hardcoded in mathias_model.m

% *************************************************************************
% Domain and discretization

% % disp('Setting up grid')
% % disp('Convention: z=0 is the BOTTOM and z=L is the TOP')

% NOTE: Convention in hydrus and Dogan and Motz is this:
% z=0 is the bottom;
% z=L is the surface.

L = 2;                  % column length [m]
Nz = 200;               % number of spatial points
z = linspace(0, L, Nz);   
dz=abs(z(2)-z(1));      % grid spacing

t_final = (3*365+10)*24*60*60;%(4*365+10)*24*60*60;  % final time [s]
Nt = 500;                % number of temporal points
t_vec = linspace(0, t_final, Nt);

% *************************************************************************

% % disp('Starting up pde solver pdepe')

m_slab = 0; % slab geometry

options=odeset('RelTol',1e-3,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7);

% Here, I define the pde, using the pdepe function in Matlab:

try

    sol = pdepe(m_slab, @(x,t,u,ux) richards_pde(x,t,u,ux,params), ...
              @(x) initial_h(x), ...
              @(zl,ul,zr,ur,t) bc_fun(zl,ul,zr,ur,t), ...
              z, t_vec,options);

    % sol is Nt x Nz x 1 (only one PDE variable u)
    h_tz=sol;

catch ME

    h_tz=zeros(length(t_vec),length(z))+100;
    sol=NaN;
    disp('pde solver blowup')

end


% % disp('done')

% *************************************************************************
% Outputs:

[Nt_temp,~]=size(h_tz);

if(Nt_temp~=Nt)
    h_tz=zeros(length(t_vec),length(z))+100;
    sol=NaN;
    disp('pde solver failed to complete')
end


theta_tz = 0*h_tz;
Kmat  = 0*h_tz;
Cmat  = 0*h_tz;



for ti = 1:Nt          
    for zi = 1:Nz
        [theta_tz(ti,zi), Cmat(ti,zi), Kmat(ti,zi)] = vg_model(h_tz(ti,zi), params);
    end
end

% OUTPUT: Re-scaling time:

t_vec=t_vec/(24*60*60); % t in days.

% *************************************************************************
% *************************************************************************

% PDEPE solves:
%
% c(du/dt)=(d/dx)(f)+s, where:
%
% c=c(x,t,u,u_x) is the capacity
% f=f(x,t,u,u_x) is the flux
% s=s(x,t,u,u_x) is the source term.

% In the Richards Equation, x is playing the role of z.   
% Similarly, u = h, u_x = dh/dz.  Notation: ux\equiv u_x.  
% % In the next subfunction, t_local is the time variable within the scope 
% of the subfunction itself, measured in seconds.

function [c,f,s] = richards_pde(x_local, t_local, u, ux, params_local) 

    % th_r = params_local.theta_r; % local value of residual moisture level.
    % th_s = params_local.theta_s; % local value of saturated moisture level.
    % z_local = x_local;      % local value of z (z\equiv x).
    h_local = u;            % local value of h (\equiv u).

    % Get theta, C, K from VG model.
    % The first output variable is theta_val.
    [theta_val, Cval, Kval] = vg_model(h_local, params_local); 

    c = Cval;             % capacity
    
    f = Kval .* (ux + 1); % flux f = K(h)*( (dh/dz) + 1); HYDRUS CONVENTION


    [f1,f2] = mathias_model(x_local,h_local,theta_val,L,params_local);     % sink term from the Mathias paper, 
                                                            % before the application of
                                                            % the EP term.                     
    s_prefactor=f1*f2;

    EP_local=get_EP(t_local,Tdata,colnumber_EP); % Obtaining EP coefficeint from data
                                                 % (depends on time).

    %     % To avoid potential blow-up of the PDE, I might need the following safety
    %     % switch:
    % 
    %     theta_upr=0.1;
    %     theta_lwr=0.09;
    % 
    %     if(theta_val>theta_upr)
    %         safety=1;
    %     elseif(theta_val>theta_lwr)
    %         safety=(theta_val-theta_lwr)/(theta_upr-theta_lwr);
    %     else
    %         safety=0;
    %     end
    % 
    %     % safety=(theta_val-th_r)/(th_s-th_r);

    s=-EP_local*s_prefactor; 

end

% *************************************************************************
% BCs for pdepe
%
% BCs have the form: p(x,t,u) + q(x,t,u)*f(x,t,u,u_x) = 0.
%
% In the code:
%
%  * zl, ul are left boundary z = 0.  Caution!  This is the BOTTOM.
%  * zr, ur are the right boundary z=L.  Caution!  This is the TOP.

function [pl,ql,pr,qr] = bc_fun(zl, ul, zr, ur, t_cur)

    % A pseudoponding condition is implemented, with q_max as the maximum
    % infiltation.

    q_max=8e-8;
    q_rainfall=get_rainfall(t_cur,Tdata,colnumber_rainfall);

    % If q_rainfall>q_max, the top level experiences excess infiltration,
    % so we switch to a Dirichlet BC (saturation at top).
    % % Otherwise, the system can accept all rainfall into the column, 
    % and we work with the standard top BC.

    if(q_rainfall>q_max) 
        pr= h_ref-ur; 
        qr=0;
        % display(strcat('ponding at day=',num2str(t_cur/(24*3600))))
    else 
        qtop=q_rainfall; 
        pr =  - qtop;
        qr = 1;
    end    
    
    % Left-hand BC (ul,rl) corresponds to z=0, which is the BOTTOM, where we 
    % enforce free drainage, hence f = K(ul).

    % [~, ~, Kbot] = vg_model(ul, params_local);
    % pl =  - Kbot;
    % ql = 1;

    % Optional - replace bottom BC with fixed prssure head.
    pl= h_ref-ul;
    ql=0;
end

% *************************************************************************

% Here, the initial condition is set.
% Example: h_init(z)=z_{wt}(0)-z.

function h0 = initial_h(x)
    % If the water table is 1 m  below the bottom of the domain, we have:
    h0 = h_ref+0*x;%-1-x;
end

% *************************************************************************
% *************************************************************************

end