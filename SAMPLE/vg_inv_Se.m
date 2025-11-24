function h = vg_inv_Se(Se, params)
%VG_INV_SE  Invert Van Genuchten effective saturation to pressure head.
%
%   h = vg_inv_Se(Se, params)
%
%   Inputs:
%       Se      - effective saturation (0 < Se <= 1)
%       params  - structure with fields:
%                   params.alpha (1/m)
%                   params.n     (dimensionless)
%
%   Output:
%       h       - pressure head (negative for unsaturated conditions)
%
%   Notes:
%       - Returns h < 0 (Richards equation convention for unsaturated soil)
%       - Se should be in (0,1]; Se=1 gives h=0 (mathematically)
%       - This function is vectorized for Se arrays

    alpha = params.alpha;
    n     = params.n;
    m     = 1 - 1/n;

    % Protect against rounding errors:
    Se = min(max(Se, eps), 1);   % Clamp Se to (0,1]

    % Compute |h|
    abs_h = ((Se.^(-1/m) - 1).^(1/n)) ./ alpha;

    % Pressure head is negative in unsaturated zone:
    h = -abs_h;
end