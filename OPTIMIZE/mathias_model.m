function [f1,f2]=mathias_model(z_hydrus,~,theta_val,L,params)


theta_r=params.theta_r;
theta_s=params.theta_s;

theta_span=theta_s-theta_r;
d_theta=theta_span/5;

theta1=theta_r+d_theta;
theta2=theta_r+2*d_theta;
theta3=theta_r+3*d_theta;
theta4=theta_r+4*d_theta;

Lr=0.25;
a=1.55;

z_mathias=L-z_hydrus;

if(z_mathias>Lr)
    f1=0;
else
    f1=(a/Lr)*(exp(-a)-exp(-a*z_mathias/Lr))/((1+a)*exp(-a)-1);
end

if(theta_val<theta1)
    f2=0;
elseif(theta_val<theta2)
    f2=(theta_val-theta1)/(theta2-theta1);
elseif(theta_val<theta3)
    f2=1;
elseif(theta_val<theta4)
    f2=1-(theta_val-theta3)/(theta4-theta3);
else
    f2=0;
end


end