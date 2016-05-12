function [smm] =smm_model_2(params, sim_number,data)

% basically, I just copy all the information I need from
% 'smm_simulation_hw3' 
%% parameters

beta=params(1);
mu_wage=params(2);
std_wage=params(3);
%gamma=params(4);
%phi=params(5);

%% observation

wage_observed=data{1};
sample_size=data{2};
age=data{3};
random_draw=data{4};
epsilon=mu_wage+std_wage.*random_draw;

%% simulation

age_sim=repmat(age, sim_number,1);
wage_simulated=beta.*age_sim+epsilon;

%% without truncation, without threshold
% threshold= gamma.*age_sim+ phi;
%home_production_index= wage_simulated<threshold;
% sim_size=sim_number-max(sum(home_production_index));

%wage_simulated(home_production_index)=NaN;
% truncated_wage=zeros(sim_size,length(age));



%% 
g=sum(wage_observed)./sample_size -sum(wage_simulated)./sim_number;
smm=g*g';

        
    








end


