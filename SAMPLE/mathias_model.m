function [f1,f2]=mathias_model(z_hydrus,h_val,L)

% Mathias et al.:
% ha=-0.05;
% hd=-4;
% hw=-150;

% From Rosetta, specialized to Johnstown Castle (SK, on 21/11/2025):

ha = -0.25; % critical pressure heads associated with anaerobiosis,
hd = -3;  % critical pressure heads associated with soilwater-limited evapotranspiration
hw = -10;  %  critical pressure head associated with plant wilting

% Caution!  Please perform the following check.
% With the given VG paramters, theta_out= 0.2014 at hw=-10, and this is
% greater than theta_r (comforatably so), so this is OK.  Otherwise, the
% solver could blow up.

Lr=0.25;
a=1.55;

z_mathias=L-z_hydrus;

if(z_mathias>Lr)
    f1=0;
else
    f1=(a/Lr)*(exp(-a)-exp(-a*z_mathias/Lr))/((1+a)*exp(-a)-1);
end

if(h_val>ha)
    f2=0;
elseif(h_val>hd)
    f2=1;
elseif(h_val>hw)
    f2=1-((h_val-hd)/(hw-hd));
else
    f2=0;
end