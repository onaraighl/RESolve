function [theta_out, C_out, K_out] = vg_model(h_val, params_vg)

% *************************************************************************
% Matlab function vg_model.m
% Computes the standard Van Genuchten model.
% Lennon Ó Náraigh
% November 2025
% 
% Inputs:
%
% * h_val - pressure head - can be positive or negative.  Positive indicates
%   saturation.
% * params_vg - VG parameters
%
% Outputs:
%
% * theta_out... soil moisture in mm^3/mm^3
% * C_out... capacity = dtheta/dh, evaluated analutically from the VG model.
% * K_out... permeability.
%
% *************************************************************************

% PARAMETERS
th_r = params_vg.theta_r;
th_s = params_vg.theta_s;
alpha = params_vg.alpha; % 1/m
n = params_vg.n;
Ks = params_vg.Ks;       % m/s
m = 1 - 1/n;
eta = params_vg.eta;
capacity_correction=params_vg.capacity_correction;

% *************************************************************************
% Return values in saturated case for h>=0:

if h_val >= 0
    theta_out = th_s;
    C_out = 1e-8; % small but nonzero to avoid zero capacity (or use eps)
    K_out = Ks;
    return;
end

% *************************************************************************
% Otherwise, work with |h| because VG defined for negative h

x = alpha * abs(h_val); % auxiliary variable, x = alpha * (-h)
Se = (1 + x^n)^(-m);    % effective saturation
theta_out = th_r + (th_s - th_r) * Se;

% Compute capacity C = dtheta/dh  (for h<0)
% Use dSe/dh=(dx/dh)(dSe/dx), hence:
% dSe/dh = m * n * alpha * x^(n-1) * (1 + x^n)^(-m-1)

if x == 0
    dSedh = m * n * alpha * 0 * (1 + 0)^(-m-1); % zero
else
    dSedh = m * n * alpha * x^(n-1) * (1 + x^n)^(-m-1);
end
C_out = (th_s - th_r) * dSedh+capacity_correction*theta_out;

% *************************************************************************
% Hydraulic conductivity using Mualem-van Genuchten
% Use eta in the prefactor as per Mathias et al.

% K_rel = (Se^eta) * [1 - (1 - Se^(1/m))^m]^2
% K = Ks * K_rel
Krel = (Se^eta) * (1 - (1 - Se^(1/m))^m)^2;
K_out = Ks * Krel;

% *************************************************************************
% Safety clamps, as recommended by none other than ChatGPT:

% if C_out < 1e-8
%     C_out = 1e-8;
% end

if K_out < 1e-14
    K_out =  1e-14;
end

% *************************************************************************

end