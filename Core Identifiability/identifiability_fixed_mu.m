function identifiability_fixed_mu
clc;close all;
clear all;
global params change e bmin
patient = load('patient39.txt');
parameters = load('p39.mat');

t = patient(:,2);                 % time vector in days
psa = patient(:,3);               % psa vector of values
androgen= patient(:,4);           % androgen vector of value

[l,~] =  size(patient);
val = 2;
cur = 1;
cycles = 3; %
maxim = 100000;
for i = 1:l
    if val~= patient(i,6)
        change(cur) = patient(i,2);
        cur = cur + 1;
        val = patient(i,6);
        cycles = cycles - 1;
        if cycles < 0
           maxim = i;
           break;
        end
    end
end
patient = patient(1:maxim,:);

params = parameters.params(1:5);

amax = max(androgen);
bmin = min(androgen);

% mu                % q2                 % d1              % gamma1          % A0    
LB(1) = 0.001;      LB(2) = 0.01;        LB(3) = 0.001;    LB(4) = 0.008;    LB(5) = amax-5;
UB(1) = 0.09;       UB(2) = bmin + 0.2;  UB(3) = 0.09;     UB(4) = 0.8;      UB(5) = amax+5;
min_mu = LB(1);
max_mu = UB(1);
x0 = [androgen(1);0.01;0.0001;psa(1);0.9*androgen(1)]; %(A,x1,x2,P,Q)

%% Find best fit
options = optimset('Algorithm','interior-point',...
    'MaxIter',5000,'MaxFunEvals',5000);
params = fmincon(@(params)Objective(params,patient(:,3)...
    ,patient(:,4),patient(:,2),x0)...
    ,params,[],[],[],[],LB',UB',[],options);
figure(17)
[t,x] = ode15s(@new_model,[0 patient(end,2)],x0,[],params); 
hold on; box on;
plot(t,x(:,1))
scatter(patient(:,2),patient(:,3));
hold off;
%maybe plot some more plot(t,x(:,1)) + data
%% Create profile-likelihood for a single parameter (mu)
% 
MU = params(1);
global fixed_mu;
fixed_mu = MU;

total_N = 100; %number of time changing mu

% Lower_factor = 0; %<1
%    Lo_mu = fixed_mu*Lower_factor;
   Lo_mu = LB(1);
% Upper_factor = 6; % > 1
%    Up_mu = Upper_factor*fixed_mu;
   Up_mu = UB(1);
step_size = (Up_mu - Lo_mu)/total_N;

Plotting_holder = zeros(2,total_N+1);
counter = 1;

global params_fixed_mu
params_fixed_mu = params;
params_fixed_mu(1) = []; %delete mu (due to fmincon).

LB(1) = [];
UB(1) = [];

    options = optimset('Algorithm','interior-point',...
    'MaxIter',5000,'MaxFunEvals',5000);

for fixed_mu = Lo_mu:step_size:Up_mu
%Find best fit alpha (while fixing beta)
params_fixed_mu = fmincon(@(params)Objective_fixed_mu(params_fixed_mu,patient(:,3)...
    ,patient(:,4),patient(:,2),x0)...
    ,params_fixed_mu,[],[],[],[],LB',UB',[],options);

Plotting_holder(:,counter) = [fixed_mu;e];
counter = counter + 1;
end
%% convert error to -LL
N = 2*length(patient(:,3)); %total data used for fitting
sigma = 1; %assume guassian additive process for data. Same variance and independent.
Plotting_holder(2,:) = N/2*log(2*pi) + N*log(sigma) + Plotting_holder(2,:)/(2*sigma^2);

%% Calculate the threshold
% psa_time_series = timeseries(patient(:,3),patient(:,2));
% psa_var = var(psa_time_series);
% and_time_series = timeseries(patient(:,4),patient(:,2));
% and_var = var(and_time_series);
% %
% length_psa = length(patient(:,3));
% length_and = length(patient(:,4));
% mean_psa = mean(patient(:,3));
% mean_and = mean(patient(:,4));
% %
% combined_mean = (length_psa*mean_psa + length_and*mean_and)/(length_psa+length_and);
% %
% first_component = (length_psa - 1)*(psa_var^2 + (mean_psa-combined_mean)^2);
% second_component = (length_and - 1)*(and_var^2 + (mean_and-combined_mean)^2);
% %
% combined_sigma_sq = (first_component + second_component)/(length_psa+length_and-1);
%
% pre_threshold = combined_sigma_sq/chi2inv(0.95,5);
pre_threshold = icdf('Chisquare',0.95,1)/2;
min_err = min(Plotting_holder(2,:));
%
threshold = pre_threshold + min_err;
%
horizon_Line = linspace(min(Plotting_holder(1,:)),max(Plotting_holder(1,:)));
threshold_L = ones(size(horizon_Line))*threshold;
%% Plot the profile-likelihood curve
figure(100)
hax=axes; 
set(0,'DefaultAxesFontSize', 16)
plot(Plotting_holder(1,:),Plotting_holder(2,:),'k','LineWidth',2)
hold on;
plot(horizon_Line,threshold_L,'--r','LineWidth',1.5);
xlim([min(Plotting_holder(1,:)) max(Plotting_holder(1,:))]);
xlabel ('\mu_m');
ylabel ('-LL');
title('Profile Likelihood for \mu_m');
line([min_mu min_mu], get(hax,'YLim'),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle',':'); %lower threshold for mu
line([max_mu max_mu], get(hax,'YLim'),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle',':'); %upper threshold for mu
aaaa=1;
end