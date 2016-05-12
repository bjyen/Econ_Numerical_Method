function m_fun=matching_fun(theta,PI,z)
global delta alpha A beta kappa mu b; 
global PI z;
% Using the surplus equation to sub out J (and a bunch of the other stuff)
% we get the following stochatic difference equation
m_fun=(kappa/(beta*A)).*theta.^(alpha)-PI*((1-mu).*(z-b)-(kappa*mu).*theta+((1-delta)*kappa/A).*theta.^(alpha));
end





