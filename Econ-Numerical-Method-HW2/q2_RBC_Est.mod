var y c k i l y_l w r z;
varexo e;

parameters beta psi delta alpha rho gamma sigma epsilon;
alpha = 0.33;
beta = 0.99;
delta = 0.023;
psi = 1.75;
rho = 0.95;
sigma = (0.007/(1-alpha));
epsilon = 10;

model;
(1/c) = beta*(1/c(+1))*(1+r(+1)-delta);
psi*c/(1-l) = w;
c+i = y;
y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha);
w = y*((epsilon-1)/epsilon)*(1-alpha)/l;
r = y*((epsilon-1)/epsilon)*alpha/k(-1);
i = k-(1-delta)*k(-1);
y_l = y/l;
z = rho*z(-1)+e;
end;

varobs y;

initval;
k = 9;
c = 0.76;
l = 0.3;
w = 2.07;
r = 0.03;
z = 0;
e = 0;
end;

// Calculate steady state and simulation
steady;
check;

shocks;
var e = sigma^2;
end;

stoch_simul(periods=2100);
// end of simulation


// Estimation
estimated_params;
alpha, beta_pdf, 0.35, 0.02;
beta, beta_pdf, 0.99, 0.002;
delta, beta_pdf, 0.025, 0.003;
psi, gamma_pdf, 1.75, 0.1;
rho, beta_pdf, 0.95, 0.05;
epsilon, gamma_pdf, 10, 0.5;
stderr e, inv_gamma_pdf, 0.01, inf;
end;
// End of estimation

estimation(datafile=q2_RBC_data,mh_replic=20000,mh_nblocks=2,mh_drop=0.45,mh_jscale=0.8,mode_compute=6);
