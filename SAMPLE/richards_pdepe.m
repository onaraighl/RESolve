function [z,t,h_tz,theta_tz,sol]=richards_pdepe()

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
% SYNTHETIC METEOROLOGICAL DATA

metData.Year=365*24*3600;
metData.RfMax=20*1e-3/(24*3600); % max daily rainfall in m/s
metData.EpMax=10*1e-3/(24*3600); % max daily ep in m/s

% *************************************************************************
% VG parameters:

% % Data for Johnstown Castle

params.theta_r = 0.077;
params.theta_s = 0.396;
params.alpha   = 0.894;     % 1/m
params.n       = 1.424;
params.Ks      = 0.195/(24*3600);  % m/s; NOTE - this is the value from Rosetta.
params.eta     = 0.5;
params.capacity_correction=9.81e-7; % From Mathias et al., under Equation (17)

Se_target=0.995;
h_ref = vg_inv_Se(Se_target, params);  % reference level, very close to saturation.

display(strcat('Reference head very close to saturation is =',num2str(h_ref),' m'))
disp('Convention: negative heads correspond to unsaturates soils')

% *************************************************************************
% Domain and discretization

disp('Setting up grid')
disp('Convention: z=0 is the BOTTOM and z=L is the TOP')

% NOTE: Convention in hydrus and Dogan and Motz is this:
% z=0 is the bottom;
% z=L is the surface.

L = 1.8;                % column length [m]
Nz = 200;               % number of spatial points
z = linspace(0, L, Nz);   
dz=abs(z(2)-z(1));      % grid spacing

t_final =365*3*24*60*60; % final time [s] ... so running simulation for 3 years.
Nt = 500;                % number of temporal points
t = linspace(0, t_final, Nt);

disp('done')

% *************************************************************************

disp('Starting up pde solver pdepe')

m_slab = 0; % slab geometry

options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7);

% Here, I define the pde, using the pdepe function in Matlab:

sol = pdepe(m_slab, @(z,t,u,ux) richards_pde(z,t,u,ux,params), ...
              @(z) initial_h(z), ...
              @(zl,ul,zr,ur,t) bc_fun(zl,ul,zr,ur,t,params,dz), ...
              z, t,options);

disp('done')

% *************************************************************************
% Outputs:

% sol is Nt x Nz x 1 (only one PDE variable u)
h_tz=sol;

theta_tz = 0*h_tz;
Kmat  = 0*h_tz;
Cmat  = 0*h_tz;

for ti = 1:Nt           
    for zi = 1:Nz
        [theta_tz(ti,zi), Cmat(ti,zi), Kmat(ti,zi)] = vg_model(h_tz(ti,zi), params);
    end
end

% OUTPUT: Flipping the coordinates and scaling the time:

% z=z(end:-1:1);  % now, z=0 is the top.
t=t/(24*60*60); % t in days.

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
% of this function, measured in seconds.

function [c,f,s] = richards_pde(x, t_local, u, ux, params_local) 

    % th_r = params_local.theta_r; % local value of residual moisture level.
    % th_s = params_local.theta_s; % local value of saturated moisture level.
    % z_local=x; % local value of z (z\equiv x).
    h_local = u; % local value of h (\equiv u).

    % Get theta, C, K from VG model.
    % The first output variable would be theta_val.
    [~, Cval, Kval] = vg_model(h_local, params_local); 

    c = Cval;             % capacity
    f = Kval .* (ux + 1); % flux f = K(h)*( (dh/dz) + 1); HYDRUS CONVENTION

    [f1,f2] = mathias_model(x,h_local,L);   % sink term from the Mathias paper, 
                                            % before the application of
                                            % the EP term.

    s_prefactor=f1*f2;

    EP_local=get_EP(t_local,metData); % Obtaining EP coefficeint from data
                                      % (depends on time).

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
%
% We enforce:
%
%  Top (z=L): specified flux q_top(t)
%             => enforce f(z=L) = q_top(t)
%             So choose: p = -q_top(t).
%
%  Bottom (z=0): free drainage => dh/dz = 0 
%             => f_bottom = K(h_bottom)*(0 + 1) = K(h_bottom)
%             So enforce f = K(ul)
%             So choose p = -K(bot).

function [pl,ql,pr,qr] = bc_fun(zl, ul, zr, ur, t_cur, params_local,dz_local)

    % Right-hand BC (ur,zr) corresponds to z=L, which is the TOP.
    % Note!  A pure ponding condition is not possible in pdepe.
    % Here, we implement the standard top BC.

    q_rainfall=get_rainfall(t_cur,metData);
    qtop=q_rainfall;

    pr =  - qtop;
    qr = 1;
    
    % % Left-hand BC (ul,rl) corresponds to z=0, which is the BOTTOM, where we 
    % % enforce free drainage, hence f = K(ul).
    % 
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
    % If the water table is 0.1 m  below the bottom of the domain, we have:
    h0 = -0.5-x;
end

% *************************************************************************
% *************************************************************************

end