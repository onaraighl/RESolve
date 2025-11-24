function [THETA]=postprocessTHETA(t,z,theta_tz)

% *************************************************************************
% Matlab function postprocessTHETA.m
% Computes the storage water level
% 
% \Theta(t)=\int_0^L \theta(t,z)dz,
%
% with values in metres.  Integration via the trapezoidal rule.
%
% Author: Lennon Ó Náraigh
% Date: November 2025
% 
% *************************************************************************
% Inputs:
%
%  t        - vector of time values
%  z        - z-coordinate
%  theta_tz - spacetime array of values of theta, in units of mm^3/mm^3.
%
% Outputs:
%
% THETA - the storage water level, as a function of t, in units of metres.
%
% *************************************************************************

dz=abs(z(2)-z(1)); % absolute value

THETA=0*t;

% Performing trapezoidal rule calculation:

for i=1:length(t)
    temp=theta_tz(i,:);
    THETA(i)=sum(temp)*dz;
end