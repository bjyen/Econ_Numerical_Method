%% Numerical Method HW1
%% Bing-Jie Yen
clear all;
clc;
%% Discrete-time version of a Mortensen-Pissarides with aggregate fluctuation
%% 1. Set parameters (global variables)
%% 2. Space discrete AR(1) : Rouwenhorst method
%% 3. Find the roots by solving stochastic difference equation  

% 1. Parameters 
global delta alpha A rho_z sigma_eps mu b beta kappa N;


delta=0.0081;             % Separation rate
alpha=0.72;               % Elasticity of matching
A=0.158;                  % Matching efficiency
rho_z=0.9895;             % Autocorrelation of weekly productivity
sigma_eps=0.0034;         % Standard Deviation for innovations
mu=0.72;                  % Bargaining weight for workappaers
b=0.4;                    % Unemployment utility
beta=0.9992;              % Weekly discount rate
kappa=0.338;              % Vacancy posting cost
N =50;                   % Number of grids
lambda = 3;


% 2. Approximate AR(1) using Rouwenhorst method for AR(1)
% (the reason why I choose to use Rouwenhorst instead of Tauchen: with
% extremely persistent processes, we oftern need many points in Tauchen and
% the perforamnce if poor)

% 2. (1) A method for choosing values for realizations and the transition
% matrix
% 2 (2) Markov Chain Approximation: 
% Purpose: Finite state Markov chain approximation of the first order
% autoregressive process AR(1)

%------hint: Consider matrix recursive form: for N=2 and N>2
%------ Normalize the matrix: at each step, normalize to be stochastic matrix by dividing by 2

% 2. (step 1): we can choose p, q to hit more rich process, but z is normal and
% we restrict grid; there is little freedom.
% Therefore, let p=q=(1+rho_z)/2 to hit autocorrelation and spacing is to
% determined to hit variance

N=50;                            % # of grids
p=0.5*(1+rho_z);                 % Initial value for p=q=(1+rho_z)/2
q=0.5*(1+rho_z);
PI=[p 1-p;1-q q];                % i=2         
 for i= 3:N
     PI= p*[PI zeros(i-1,1);zeros(1,i-1) 0]+ (1-p)*[ zeros(i-1,1) PI; 0 zeros(1,i-1)]+(1-q)*[zeros(1,i-1) 0;PI zeros(i-1,1)] + q*[zeros(1,i-1) 0 ; zeros(i-1,1) PI];
     % At each step, normalize to be a stochastic matrix divide by 2 (NOTE:
     % NOT include first and last row
     PI(2:end-1,:)=PI(2:end-1,:)/2;
 end
 
% 2(step2): itwill closely mimic the underlying continuous valued autoregressive
 % process. z_(t+1)=rho_z*z_t+epsilon ~N(0, sigma_eps^2) Therefore,
 
sigma_z = sigma_eps/sqrt(1-rho_z^2);
z_low=-lambda*sigma_z;
z_up = lambda *sigma_z;
z = linspace(z_low,z_up,N)'; 
step= (-2*z_low)/(N-1);
for i=1:N
    z(i)= z_low + (i-1)*step
end

% 2.(2) Solve N non-linear equations: 
% Matching function= m(u,v)= A*(u^alpha)*(v^(1-alpha))
% p_theta= m/v A.*theta.^(1-alpha); job finding rates(offers arrive at
% rate)
% q_theta=theta*p_theta     Applications arrive, the rate at which vacancies are filled         
% where market tightness: theta=v/u
global PI z;

theta0=zeros(N,1); % Guess of initial values
theta=fsolve(@matching_fun,theta0); % Solution theta of the non-linear system
p_theta=A.*theta.^(1-alpha);
q_theta=A.*theta.^(-alpha);

 % If p_theta(i)>1, it means the market tightness is the same as previous
 % period, theta(i-1)=thta(i)
 % If q_theta(i)>1, it means the market tightness in the next period is the
 % same as the current period, theta (i)= theta(i+1)
for i=2:N;
    if p_theta(i)>1
        theta(i)=theta(i-1)
    end
end
for i=1:N-1
    if q_theta(i)>1
        theta(i)=theta(i+1)
    end
end
     

subplot(3,1,1)
plot(z,theta)
xlabel('productivity shock z');
ylabel('theta');


subplot(3,1,2)
plot(z,q_theta)
xlabel('productivity shock z');
ylabel('q_theta');


subplot(3,1,3)
plot(z,p_theta)
xlabel('productivity shock z');
ylabel('p_theta');

% 3. Simulate some realizations of the productivity and compute
% time-series of the endogenous variables
% 3(1) Simulating from the Markov Chain
% Step1: calculate the eigenvectors and correpsonding eigenvalues of PI
% Step 2: Normalize eigenvector
% Step 3: set a 'seed' so the results would be the same in each round
% Step 4: Here I use the following method to initialize the starting state
% (1) Initialize from the unconditional distribution
% (2) Use PI to get the unconditional distribution: pi
% (3) Draw a unifrom random variable x [0,1]
% (4) set z_0=z_i such that z<= sum pi

[V,D]=eig(PI);
pi=V(:,1)/sum(V(:,1));
pi_cumulated=cumsum(pi);

T=1500;
s=rng;  %set a 'seed' Save the current state of the random number generator and create a T-by-1 vector of random numbers.
x=rand(T,1);
sim= nan(T,1);
sim(1)= find((cumsum(pi)>=x(1)),1,'first');
        for j=1:T-1
            sim(j+1)=find((cumsum(PI(sim(j),:))>=x(j+1)),1,'first');
        end
        
theta_sim=theta(sim);
p_sim=A.*theta_sim.^(1-alpha);
q_sim=A.*theta_sim.^(-alpha);

subplot(3,1,1)
plot(1:T,theta_sim)
xlabel('Number of Periods');
ylabel('theta_sim');

subplot(3,1,2)
plot(1:T, q_sim)
xlabel('Number of Periods');
ylabel('q_theta_sim');

subplot(3,1,3)
plot(1:T,p_sim)
xlabel('Number of Periods');
ylabel('p_theta_sim');
            


clear all;
%% Hagedorn and Manovskii (2008) solved the same thing with the following calibrartion 
%mu=0.05
% b=0.95

% 1. Parameters 
global delta alpha A rho_z sigma_eps mu b beta kappa N;


delta=0.0081;             % Separation rate
alpha=0.72;               % Elasticity of matching
A=0.158;                  % Matching efficiency
rho_z=0.9895;             % Autocorrelation of weekly productivity
sigma_eps=0.0034;         % Standard Deviation for innovations
mu=0.05;                  % Bargaining weight for workappaers
b=0.95;                    % Unemployment utility
beta=0.9992;              % Weekly discount rate
kappa=0.338;              % Vacancy posting cost
N =50;                   % Number of grids
lambda = 3;


% 2. Approximate AR(1) using Rouwenhorst method for AR(1)
% (the reason why I choose to use Rouwenhorst instead of Tauchen: with
% extremely persistent processes, we oftern need many points in Tauchen and
% the perforamnce if poor)

% 2. (1) A method for choosing values for realizations and the transition
% matrix
% 2 (2) Markov Chain Approximation: 
% Purpose: Finite state Markov chain approximation of the first order
% autoregressive process AR(1)

%------hint: Consider matrix recursive form: for N=2 and N>2
%------ Normalize the matrix: at each step, normalize to be stochastic matrix by dividing by 2

% 2. (step 1): we can choose p, q to hit more rich process, but z is normal and
% we restrict grid; there is little freedom.
% Therefore, let p=q=(1+rho_z)/2 to hit autocorrelation and spacing is to
% determined to hit variance

N=50;                            % # of grids
p=0.5*(1+rho_z);                 % Initial value for p=q=(1+rho_z)/2
q=0.5*(1+rho_z);
PI=[p 1-p;1-q q];                % i=2         
 for i= 3:N
     PI= p*[PI zeros(i-1,1);zeros(1,i-1) 0]+ (1-p)*[ zeros(i-1,1) PI; 0 zeros(1,i-1)]+(1-q)*[zeros(1,i-1) 0;PI zeros(i-1,1)] + q*[zeros(1,i-1) 0 ; zeros(i-1,1) PI];
     % At each step, normalize to be a stochastic matrix divide by 2 (NOTE:
     % NOT include first and last row
     PI(2:end-1,:)=PI(2:end-1,:)/2;
 end
 
% 2(step2): itwill closely mimic the underlying continuous valued autoregressive
 % process. z_(t+1)=rho_z*z_t+epsilon ~N(0, sigma_eps^2) Therefore,
 
sigma_z = sigma_eps/sqrt(1-rho_z^2);
z_low=-lambda*sigma_z;
z_up = lambda *sigma_z;
z = linspace(z_low,z_up,N)'; 
step= (-2*z_low)/(N-1);
for i=1:N
    z(i)= z_low + (i-1)*step
end

% 2.(2) Solve N non-linear equations: 
% Matching function= m(u,v)= A*(u^alpha)*(v^(1-alpha))
% p_theta= m/v A.*theta.^(1-alpha); job finding rates(offers arrive at
% rate)
% q_theta=theta*p_theta     Applications arrive, the rate at which vacancies are filled         
% where market tightness: theta=v/u
global PI z;

theta0=zeros(N,1); % Guess of initial values
theta=fsolve(@matching_fun,theta0); % Solution theta of the non-linear system
p_theta=A.*theta.^(1-alpha);
q_theta=A.*theta.^(-alpha);

 % If p_theta(i)>1, it means the market tightness is the same as previous
 % period, theta(i-1)=thta(i)
 % If q_theta(i)>1, it means the market tightness in the next period is the
 % same as the current period, theta (i)= theta(i+1)
for i=2:N;
    if p_theta(i)>1
        theta(i)=theta(i-1)
    end
end
for i=1:N-1
    if q_theta(i)>1
        theta(i)=theta(i+1)
    end
end
     

subplot(3,1,1)
plot(z,theta)
xlabel('productivity shock z');
ylabel('theta');


subplot(3,1,2)
plot(z,q_theta)
xlabel('productivity shock z');
ylabel('q_theta');


subplot(3,1,3)
plot(z,p_theta)
xlabel('productivity shock z');
ylabel('p_theta');

% 3. Simulate some realizations of the productivity and compute
% time-series of the endogenous variables
% 3(1) Simulating from the Markov Chain
% Step1: calculate the eigenvectors and correpsonding eigenvalues of PI
% Step 2: Normalize eigenvector
% Step 3: set a 'seed' so the results would be the same in each round
% Step 4: Here I use the following method to initialize the starting state
% (1) Initialize from the unconditional distribution
% (2) Use PI to get the unconditional distribution: pi
% (3) Draw a unifrom random variable x [0,1]
% (4) set z_0=z_i such that z<= sum pi

[V,D]=eig(PI);
pi=V(:,1)/sum(V(:,1));
pi_cumulated=cumsum(pi);

T=2500;
s=rng;  %set a 'seed' Save the current state of the random number generator and create a T-by-1 vector of random numbers.
x=rand(T,1);
sim= nan(T,1);
sim(1)= find((cumsum(pi)>=x(1)),1,'first');
        for j=1:T-1
            sim(j+1)=find((cumsum(PI(sim(j),:))>=x(j+1)),1,'first');
        end
        
theta_sim=theta(sim);
p_sim=A.*theta_sim.^(1-alpha);
q_sim=A.*theta_sim.^(-alpha);

subplot(3,1,1)
plot(1:T,theta_sim)
xlabel('Number of Periods');
ylabel('theta_sim');

subplot(3,1,2)
plot(1:T, q_sim)
xlabel('Number of Periods');
ylabel('q_theta_sim');

subplot(3,1,3)
plot(1:T,p_sim)
xlabel('Number of Periods');
ylabel('p_theta_sim');
            







