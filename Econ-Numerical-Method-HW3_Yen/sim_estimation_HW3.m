%% Numerical Methods  HW3
%% Bing-Jie Yen
%% useful data link=>NIPA: http://www.bea.gov/iTable/index_nipa.cfm    ; FRED: research.stlouisfed.org/fred2/

%% **********Question 1***********%%
%% Parameter setting : the following is my imaginary number, and I want to see whether my estimation method would lead me to these numbers
beta     = 30;
mu_wage  = 0;
std_wage = 300;
sigma_wage= std_wage^2;
gamma     = 20;
phi       = 5;
%% Construct the observed data: 
%% step1:  Assume the working age is from 16 to 65, N=200
%% step2:  set up the labor suply function => wage=beta*t+epsilon; epsilon=mu_wage+ std_wage*rand(sample_size, length(age))

age =[16:65];
sample_size=200;
epsilon=mu_wage+std_wage*randn(sample_size,length(age));

age_sample=repmat(age,200,1); % i.e. map row vector (age) for 200 rows and one column , so dim 200x50
wage=beta.* age_sample +epsilon;

%% Truncation: set up the threshold wage between work and home production, if wage< truncation wage, choose home producation
% note: home_production_index is a indicator variable

threshold= gamma*age_sample+phi;
home_production_index = wage <threshold;

truncate_sample_size=sample_size-max(sum(home_production_index)); 
wage(home_production_index)=NaN;     %I assume home production is not categorized in labor income

truncate_wage= zeros(truncate_sample_size,length(age)); 
            % eg. BJ is 21yr old, the truncate_wage(BJ,1): the wage BJ required as she has one year working experience
            % eg. of course for BJ, the maximum workinf period is 65-21=
            % 44, which means truuncate_wage(BJ, 45:65)= []

for i=1:length(age)
    temp=wage(:,i);
    temp(isnan(temp))=[]; % if there is no wage info, return missing
    truncate_wage(:,i)= temp(1:truncate_sample_size);
end

%% SMM¡@estimate%%
%% here, we start to estimate by starting with guess values for beta, mu_wage, std_wage, gamma and phi

beta_guess      =20; % observed 30
mu_wage_guess   =5;  %0
std_wage_guess  = 100;  %300
gamma_guess     =10; %20
phi_guess       =1 ;  % 5

tau=100;
sim_number=tau.*sample_size; % every sample repeats 100 times

%% Remember that the random draw for every parmeters coming should be from the same random draw

random_draw=randn(sim_number, length(age));

%% Setting up:   data(truncate_wage, truncate_sample_size, age, randon draw)
%% guess=(5 parameters); 
%% options=optimset('Display','FunValCheck','MaxFunEvals','MaxIter','OutputFcn','PlotFcns','TolFun','TolX'); 
%% reference: http://www.mathworks.com/help/matlab/math/setting-options.html
%% [x, fval,exitflag, output]= fminsearch(function, initial value,
%% options, sim_number, data)

data={truncate_wage, truncate_sample_size, age, random_draw};
guess=[beta_guess, mu_wage_guess, std_wage_guess,gamma_guess,phi_guess];
options=optimset('MaxFunEvals',10^5,'MaxIter',10^3,'Display','iter','TolFun',1e-15,'TolX',1e-15);
[params_1,fval_1,exitflag_1]=fminsearch('smm_model_1',guess, options,sim_number,data);

params_1 
%params_1 =

 %  29.2790    1.3238  253.4041   24.0425    1.5577
 
%% Without Truncation
[params_2, fval_2, exitflag_2]=fminsearch('smm_model_2',guess, options,sim_number,data);

%params_2 
%   27.1290  181.6563 -528.2053 -411.6592 -836.5874

%% ***** Question 2: A model of separation*****%%
%Now, the same workers will face some search frictions. Job probability
%  f(t)=phi*t+lambda
%  propose and recover parameters phi, and lambda
% Background:I set up the working age is from 31 to 40 (10 periods), so the
% match will automatically terminate after 10 period. i.e get a job from 31
% and retire at 40


age_separation=[31:40];
age_sample=repmat(age_separation,200,1);
mu_productivity= 10;
std_productivity=5;
productivity_1b=50;  % the lower bound of productivity, i.e. if productivity is lower than this, get fired

productivity= mu_productivity.* age_sample + std_productivity.* age_sample.*randn(sample_size,length(age_separation));
test_index= productivity< productivity_1b;

% Note that if a worker gets separated at period t, then he would not be
% separated again for period n>t
% I construct a matrix to specify the exact period this worker gets
% separated

separation_index = zeros(sample_size,1)-1;
for i= 1:sample_size
    temp=productivity(i,:);
    if isempty(find(temp < productivity_1b,1))
        separation_index(i)=0;
    else
        separation_index(i)=find(temp<productivity_1b,1);
    end
end

new_separation= productivity< productivity_1b;
for i=1:sample_size
    if separation_index(i)==0 ||separation_index(i)==10
        continue
    else
        modify_index=separation_index(i)+ 1:length(age_separation);
        
        new_separation(i,modify_index)=0;
    end
end

    separation_rate=sum(new_separation)/sample_size;
    
%% SMM : a model of separation

mu_guess=9;
std_guess=4;
productivity_lb_guess=45;

random_draw=randn(sim_number,length(age_separation));
data ={separation_rate,sample_size,age_separation,random_draw};

guess= [mu_guess,std_guess,productivity_lb_guess];
[params_3,fval_3,exitflag_3]=fminsearch('smm_model_3',guess,options,sim_number,data);

%params_3 =

 %   3.8147    1.0136   61.6970

sim_number=sim_number*10;

random_draw=randn(sim_number,length(age_separation));
data={separation_rate,sample_size, age_separation,random_draw}
[params_4,fval_4,exitflag_4]=fminsearch('smm_model_3',guess,options,sim_number,data);

%params_4 =

%    8.5884    4.2336   46.7344
 

        












    








%% truncate_wage=
%% SMM estimate
%% Same random draw for different parameters
%% SMM estimate, ignore selection

%%
