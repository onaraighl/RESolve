function []=postprocessPlot()

% load data from file
temp2=load('data.mat');
t=temp2.temp.t; % time in days
z=temp2.temp.z; % z coordinate
h_tz=temp2.temp.h_tz; % spacetime array of pressure head, in SI units.
theta_tz=temp2.temp.theta_tz; % spacetime array of theta values, in units of volume/volume.

[THETA]=postprocessTHETA(t,z,theta_tz);

t_years=t/365; % time in years

tt=t_years+1;

plot(tt,THETA/1e-3,'linewidth',2,'color','blue') % THETA in mm (not metres)
grid on
xlim([2 4]) % start plot in 2.
ylabel('\Theta (mm)')
xlabel('Year')
% ylim([950 1420])
title('PDEPE simulation, SAMPLE')

