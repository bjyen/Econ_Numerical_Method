function [smm] = smm_model_3(params,sim_number,data)


%parameters
mu              =params(1);
std             =params(2);
productivity_lb =params(3);

%data
separation_rate =data{1};
sample_size= data{2};
age_separation=data{3};
random_draw=data{4};

%simulation

age_sim=repmat(age_separation,sim_number,1);
productivity= mu.*age_sim+std.*age_sim.*random_draw;
test_index= productivity< productivity_lb;

%% separation index

separation_index=zeros(sim_number,1)-1;

for i =1:sim_number;
    temp=productivity(i,:);
    if isempty(find(temp < productivity_lb,1))
        separation_index(i)=0;
    else
        separation_index(i)=find(temp<productivity_lb,1);
    end
end
        
%% new separation index

new_separation_index =productivity<productivity_lb;

for i=1: sim_number
    if separation_index(i)==0 || separation_index(i)==10;
        continue
    else
        modify_index=separation_index(i)+1:length(age_separation);
        new_separation_index(i,modify_index)=0;
    end
end

separation_rate_sim=sum(new_separation_index)./sim_number;
g= separation_rate - separation_rate_sim;
smm=g*g';
        







