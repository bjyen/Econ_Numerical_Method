function [smm] =smm_model_1(params, sim_number,data)

% basically, I just copy all the information I need from
% 'smm_simulation_hw3' 
%% parameters

beta=params(1);
mu_wage=params(2);
std_wage=params(3);
gamma=params(4);
phi=params(5);

%% observation

wage_observed=data{1};
sample_size=data{2};
age=data{3};
random_draw=data{4};
epsilon=mu_wage+std_wage.*random_draw;

%% simulation

age_sim=repmat(age, sim_number,1);
wage_simulated=beta.*age_sim+epsilon;

threshold= gamma.*age_sim+ phi;
home_production_index= wage_simulated<threshold;
sim_size=sim_number-max(sum(home_production_index));

wage_simulated(home_production_index)=NaN;
truncated_wage=zeros(sim_size,length(age));

if sim_size==0
    truncated_wage=truncated_wage+inf;
else
    for i=1:length(age)
        temp=wage_simulated(:,i);
        temp(isnan(temp))=[];
        truncated_wage(:,i)=temp(1:sim_size);
    end
end

%% 
g=sum(wage_observed)./sample_size -sum(truncated_wage)./sim_size;
smm=g*g';

        
    








end

