function forecasting_62_trial10
%{
    This function saves the fitting and forecasting errors for all models
    and saves the plots forecasting plots. There is one plot per (usable) patient
    and the plots include all three model forecasting.

    To save errors, creat folders 'two_pop', 'three_pop', and 'portz'.
    TO save plots,  create folder 'Plots'
    Choose patients in ~line 79
%}

%%%%%%%%%%%%%%%%%%%%%%%% Model Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose whether you want to fit all patients (q=1) or not (q=0)
q =2;
if q == 0
    LIST=[1,6,12,14,15,16,17,19,24,25,28,29,31,32,36,37,39,40,42,44,51,52,54,55,58,60,61,62,63,...
        64,66,75,77,78,79,83,84,85,87,88,91,93,94,95,96,97,99,100,101,102,104,105,106,107,108,109]; 
%     LIST = [1,6,12,14:17,19,24:25,28,29,31,32,36:37,39:40,42,44,51:52,54,55,58,60:64,66,75 ...
%         ,77:79,83:85,87,88,91,93:97,99:102,104:109]; % patient numbers 
elseif q == 1
%      LIST=[1,6,14,15,17,19,28,29,32,36,37,39,51,52,54,55,60,62,63,66,75,77,79,93,99,100,105]; 
     LIST=[1,6,14,15,17,19,28,29,32,36,37,39,51,52,54,55,60,62,63,66,75,77,79,93,100,105]; 
    %1,6,14,15,17,19,25,28,29,32,36,37,39,51,52,54,55,60,62,63,66,75,77,79,93,99,100,105,109
%     LIST = 109;
elseif q == 2
%     LIST=[1,15];
%     LIST =[1,6,14,15,17,19,25,28,29,32,36,37,39,51,52,54,55,60,62,63,66,75,77,79,93,100]; 
% LIST = [32,36,37,39];% 17 28 29 75 79];
    LIST = [15];
else
    warning('Unexpected q value, choose a different one')
    return;
end

% Control Parameters 
nFitting = 3;                       % Number of cycles used for fitting data
nForecast = 2;                      % Number of cycles of data to forecast
total_n = nFitting + 1 + nForecast; % Total number of treatment cycles 

fittingErrorPsa_two_pop = zeros(length(LIST));  
forecastErrorPsa_two_pop = zeros(length(LIST)); 
fittingErrorAnd_two_pop = zeros(length(LIST)); 
forecastErrorAnd_two_pop = zeros(length(LIST)); 

fittingErrorPsa_new = zeros(length(LIST));  
forecastErrorPsa_new = zeros(length(LIST)); 
fittingErrorAnd_new = zeros(length(LIST)); 
forecastErrorAnd_new = zeros(length(LIST)); 

fpath = 'Plots/'; % Plot destination

% Run through models
for yy = [1,2]
    if yy == 1
        model = 'two_pop';
    elseif yy == 2
        model = 'newTwo_pop';
    end

    patientIndex = 1; % initialize index to store patients who are forecasted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Run through all patients   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = LIST
    
    switch model    % Load parameters    
        case 'two_pop'            
            cd parameters_two_pop
            file_to_load = strcat('p',num2str(kk),'.mat');
            load(file_to_load);
%             load('p1.mat')
            cd ..
            par = params';
        case'newTwo_pop'
            cd parameters_newTwo_pop
            file_to_load = strcat('p',num2str(kk),'.mat');
            load(file_to_load);
%             load('p1.mat')
            cd ..
            par = params';
        otherwise 
            warning('Unexpected model, choose a different one')
            return;
    end
    
%     patientIndex = 1; % initialize index to store patients who are forecasted
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Run through all patients   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for kk = LIST
        %change = zeros(1, total_n);                      % Initialize change vector that saves time array index of when treatment changes
        % Loads Patient data for each patient in list
        patient   = strcat('patient',num2str(kk));       % Patient with corresp number
        file      = strcat('Data_62/',patient,'.txt');      % Complete name of file patient#.txt
        var       = load(file);                          % holds variable just loaded
        patient   = var;
        t         = patient(:,2);                        %time
        psa       = patient(:,3);                        %psa data
        androgen  = patient(:,4);                        %androgen data
        treatment = patient(:,6);                        %1 is on 0 is off
        index = find(LIST==kk);
        
        
        % Assign values to change vector
%         jj = 1;
%         change(1) = 1;                        % Treatment starts at t = 0
        jj = 0;
        for a = 1:length(treatment)         
            if a==1
                jj = jj + 1;
                change(jj) = a;           % Stores time in change vector 
            elseif treatment(a) - treatment(a-1) ~= 0
                jj = jj + 1;
                change(jj) = a;           % Stores time in change vector
            elseif treatment(a) ~= mod(jj,2) % When treatment change occurs
                jj = jj + 1;
                change(jj) = a;           % Stores time in change vector
            end
        end
        change(jj+1) = length(treatment);     % Last day of treatment
        
        % Determine if there is only one data point for each treatment cycle. b=1 if more than one, b=0 if not
        b = 1;                  % initialize b value
        for i = 1:length(change)-1
            a = change(i+1)-change(i);       % a(i) tells you lngth of the ith treatment
            if a == 1
                b = 0;
                break;
            end
        end   
        % Do forecasting if there are enough cycles for forecasting and if each cycle has more than one data point
        if (length(change) >= total_n) && b           
            patientNumber(patientIndex) = kk;       % Saves the patients who are forcasted
            patientIndex = patientIndex + 1;
            %Part of plotting
            figure(kk);
            box on;
            title(strcat('patient',num2str(kk)), 'fontsize',26);
            hold on;
            
            %get errors
            switch model
                case 'two_pop'
%                     params = par(:,index)';
                    params = par';
%                     x0 = [androgen(1);par(20,index);par(21,index);psa(1)]; %collects init cond in a column vector
                    x0 = [androgen(1);params(18);params(19);psa(1)];
                    [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    errp = errors(1);
                    erra = errors(2);
                    
                    [errorsf] = future(params,psa,androgen,t,change,x0_new,nFitting,nForecast,model);     % get forecasting errors and plot the forecasting
                    errpf = errorsf(1);
                    erraf = errorsf(2);

                    fittingErrorPsa_two_pop(index) = errp;
                    forecastErrorPsa_two_pop(index) = errpf;
                    fittingErrorAnd_two_pop(index) = erra;
                    forecastErrorAnd_two_pop(index) = erraf;       
    
                    case 'newTwo_pop'
%                     params = par(:,index)';
                    params = par';
                    x0 = [androgen(1);params(11);params(12);psa(1);0.9*androgen(1)];
                            [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    
                    [errors,x0_new] =  objective(params,psa,androgen,t,change,x0,nFitting,model);     % get new initial conditions, fitting errors, and plot fitting
                    errp = errors(1);
                    erra = errors(2);
                    
                    [errorsf] = future(params,psa,androgen,t,change,x0_new,nFitting,nForecast,model);     % get forecasting errors and plot the forecasting
                    errpf = errorsf(1);
                    erraf = errorsf(2);     % get forecasting errors and plot the forecasting                    

                    fittingErrorPsa_new(index) = errp;
                    forecastErrorPsa_new(index) = errpf;
                    fittingErrorAnd_new(index) = erra;
                    forecastErrorAnd_new(index) = erraf; 
            end
            %Part of plotting
            scatter(t(1:change(total_n)),psa(1:change(total_n)),'k','linewidth',3);        % plot trial data
            hold on;
            line([t(change(nFitting+1)),t(change(nFitting+1))],ylim,'linewidth',2);         % plot line that separates the fitting and forecasting sections of the plot
            xlabel('Time (Days)', 'fontsize', 26);
            ylabel('PSA (\mug/L)', 'fontsize', 26);
            p1 = plot(nan,nan,'m');
            p2 = plot(nan,nan,'g');
            p3 = plot(nan,nan,'b');
            legend([p1 p2 p3],'BK model','New model','Location','NorthWest')
            hold off;
            fileNamePSA = strcat('Patient',num2str(kk),'PSAForecast');   % Create string for androgen plot file name
            savefig(fullfile(fpath, fileNamePSA));
            saveas(gcf, fullfile(fpath, fileNamePSA), 'jpeg');
            
            
        end
            
    end
    
%     for w = patientIndex
%         figure(w)
%         fileNameAndrogen = strcat('Patient',num2str(kk),'Androgen');   % Create string for androgen plot file name
%         saveas(gcf, fullfile(fpath, fileNameAndrogen), 'jpeg');
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Error Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch model
    
        case 'two_pop'
            save('two_pop/fittingErrorPsa_two_pop.mat','fittingErrorPsa_two_pop') 
            save('two_pop/fittingErrorAnd_two_pop.mat','fittingErrorAnd_two_pop')
            save('two_pop/forecastErrorPsa_two_pop.mat','forecastErrorPsa_two_pop')
            save('two_pop/forecastErrorAnd_two_pop.mat','forecastErrorAnd_two_pop')
        case 'newTwo_pop'
            save('new/fittingErrorPsa_new.mat','fittingErrorPsa_new') 
            save('new/fittingErrorAnd_new.mat','fittingErrorAnd_new')
            save('new/forecastErrorPsa_new.mat','forecastErrorPsa_new')
            save('new/forecastErrorAnd_new.mat','forecastErrorAnd_new')
    end
end
end

function [errors,x0_new] = objective(params,psadata,and_data,tdata,change,x0,n,model) 
%{
 This function takes the parameters from the fitting and the initial
 conditions from the trial data.

 This function plots the fitting model data and returns the fitting 
 errors and the values of the dependent variables at the end of the
 fitting cycle.
%}

psadata = psadata(change(1):change(n+1));    
and_data = and_data(change(1):change(n+1));
[Y1run] = run_model(params,tdata,change,x0,1,n,model,and_data);

switch model
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        errors = [errp, erra];
        
        %Update initial conditions
        X10_new = Y1run(end,2);
        X20_new = Y1run(end,3);
        PSA0_new = Y1run(end,4);
        Q0_new = Y1run(end,1);
        
        x0_new = [Q0_new, X10_new, X20_new, PSA0_new];
        
        plot(tdata(1:length(PSA)),PSA,'m','linewidth',2);
        hold on;

    case 'newTwo_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);
        errp = sum((PSA-psadata).^2/length(PSA));
        erra = sum((AND-and_data).^2/length(AND));
        errors = [errp, erra];
        
        %Update initial conditions
        X10_new = Y1run(end,2);
        X20_new = Y1run(end,3);
        PSA0_new = Y1run(end,4);
        Q0_new = Y1run(end,5);
        AND_new = Y1run(end,1);
        x0_new = [AND_new, X10_new, X20_new, PSA0_new, Q0_new];
        
        plot(tdata(1:length(PSA)),PSA,'g','linewidth',2);
        hold on;
end


    
end

function [errors]= future(params,psadata,and_data,tdata,change,x0,nFitting,nForecast,model) 
%{
 This function takes the parameters from the fitting and the new initial conditions
 from where the fitting model left off to run the model again.

 It adds a forecasting plot and returns the forecasting errors.
%}
psadata = psadata(change(nFitting+1):change(nFitting+1+nForecast));                % retrieve psa data from after fitting cycles
androgen = and_data;    % Androgen used for portz
and_data = and_data(change(nFitting+1):change(nFitting+1+nForecast));       %androgen used for errors                          % retrieve androgen data from after the fitting cycles
[Y1run] = run_model(params,tdata,change,x0,nFitting+1,nFitting+nForecast,model,androgen);             % run model with the initial conditions set to where the fitting left off
switch model
    case 'two_pop'
        PSA = Y1run(:,4);
        AND = Y1run(:,1);        
        errp = sum((PSA-psadata).^2/length(PSA));                        % forecasting psa error
        erra = sum((AND-and_data).^2/length(AND));                       % forecasting androgen error
        errors = [errp, erra];
        
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'m','linewidth',2);            % plot the forecasting model data
        hold on;
        
     case 'newTwo_pop'
     PSA = Y1run(:,4);
        AND = Y1run(:,1);        
        errp = sum((PSA-psadata).^2/length(PSA));                        % forecasting psa error
        erra = sum((AND-and_data).^2/length(AND));                       % forecasting androgen error
        errors = [errp, erra];
        
        plot(tdata(change(nFitting+1):change(nFitting+1+nForecast)),PSA,'g','linewidth',2);            % plot the forecasting model data
        hold on;   
end

end

%%%%%%%%%%%%%%%%  Runs the Model %%%%%%%%%%%%

function [y] = run_model(params,tdata,change,x0,initial_n, final_n, model,androgen)   
y = [];                                      % Initial vector for solution
tint = tdata(change(final_n+1));             % Time interval to run solution 
for k = initial_n:final_n                    % K number of cycles of treatment
    u = 1 - mod(k,2);    
    switch model             
        case 'two_pop'
            
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            
           if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0,options);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)]; %#ok<AGROW>
                end
            end
            
         case 'newTwo_pop'
                options = odeset('RelTol',1e-12,'AbsTol',1e-14);
            if tint(end) >= tdata(change(k))
                [~,Yrun]=ode15s(@(t,x) newTwo_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0,options);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < final_n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                elseif k == final_n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                end
            end
%            
    end   
end
end

%% Model ODE functions

function dxdt = two_pop(~,x,p)
% collect parameter values to pass to ODE function
um = p(1);        q1= p(2);        q2= p(3);      c1= p(4);    K1= p(5);
b= p(6);         sigma1= p(7);  epsilon= p(8);   d1= p(9);  d2= p(10);
R1= p(11);        R2= p(12);       gamma1= p(13); gamma2= p(14);    dd1 = p(15);
dd2 = p(16); Qm = p(17); u= p(20);
%separates solutions
Q = x(1); X = x(2); Y = x(3); P = x(4);
%%%%%%%%%%%%%%%%%%%  Parameters for ODE System %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q1>Q
    ux = 0;
else
    ux = um*(1 - q1/Q);
end
if q2 > Q
    uy = 0;
else
    uy = um*(1 - q2/Q);
end
Dx = d1*R1/(Q+R1);    Dy = d2*R2/(Q+R2);
mxy = c1*K1/(Q + K1); 
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = (ux - Dx  - dd1*X)*X -  mxy*X;
dY = (uy - Dy  - dd2*Y)*Y + mxy*X;
dA = (gamma1*u + gamma2)*(Qm -Q) - (ux*Q*X + uy*Q*Y)/(X+Y);
dP = b*Q + sigma1*(Y*Q + X*Q) - epsilon*P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY;dP]; %Puts the ode in a column vector
end


function dxdt = newTwo_pop(~,x,p) 
% global d2 R2 gamma2 c K
% collect parameter values to pass to ODE function
u1 = p(1); u2 = p(2); q1= p(3); q2= p(4);
b = p(5); d1 = p(6); R1 = p(7); gamma1 = p(8);
sigma = p(9); epsilon = p(10);
A0 = p(13);  
dd1 = 5; dd2 = 5;
u= p(17);
%
d2 = p(14); R2 = p(15); gamma2 = p(16);
c = 0.00015; K = 1;
%separates solutions 
Q = x(5); X = x(2); Y = x(3); P = x(4); A = x(1);
%%%%%%%%%%%%%%%%%%%  Parameters for ODE System %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q1>Q
    ux = 0;
else
    ux = u1*(1 - abs(q1)/Q);
end
if q2 > Q
    uy = 0;
else
    uy = u2*(1 - abs(q2)/Q);
end
Dx = d1*R1/(Q+R1);    Dy = d2*R2/(Q+R2);
mxy = c*K/(Q + K);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = ux*X - Dx*X - dd1*X*X - mxy*X;
dY = uy*Y - Dy*Y - dd2*Y*Y + mxy*X;
dA = gamma1*(A0 - A) - gamma1*A0*(1-u);
dQ = gamma2*(A - Q) - (ux*Q*X + uy*Q*Y)/(X+Y);
dP = b*Q + sigma*(Y*Q + X*Q) - epsilon*P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY; dP; dQ]; %Puts the ode in a column vector
end