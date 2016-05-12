//  The stochastic case
// Refer User Guide Chapter 3
// You can find the complete code : 3.9.1 The stochastic model

var y c k i l y_l w r z;
varexo e;  // with shocks

parameters beta psi delta alpha rho sigma epsilon;
alpha = 0.33;
beta = 0.99;
delta = 0.023;
psi = 1.75;
rho = 0.95;
sigma = (0.007/(1-alpha));
epsilon = 10;


/*3.5.1 Model in Dynare notation*/
/*Just in case you need a hint or two to recognize these equations, here’s
a brief description: the ﬁrst equation is the Euler equation in consumption.
The second the labor supply function. The third the accounting identity. The
fourth is the production function. The ﬁfth and sixth are the marginal cost
equal to markup equations. The seventh is the investment equality. The
eighth an identity that may be useful and the last the equation of motion of
technology.*/

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

/* 3.6.1 Stochastic models and steady states*/

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

/* 3.7.3 Stochastic models*/
shocks;
var e = sigma^2;
end;

stoch_simul(periods=2100);
// end of simulation

/* User Guide chapter 5
// Estimation 


estimated params;
alpha, beta pdf, 0.35, 0.02;
beta, beta pdf, 0.99, 0.002;
delta, beta pdf, 0.025, 0.003;
psi, gamma pdf, 1.75, 0.1;
rho, beta pdf, 0.95, 0.05;
epsilon, gamma pdf, 10, 0.5;
stderr e, inv gamma pdf, 0.01, inf;
end;
// End of estimation

estimation(datafile=simuldataRBC,nobs=200,first obs=500,
mh replic=2000,mh nblocks=2,mh drop=0.45,mh jscale=0.8,
mode compute=4);


