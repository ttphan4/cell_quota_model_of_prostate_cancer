function [out] =  trial10
%{
This files takes the patients numbers from variable LIST and performs parameter
estimation using fmincon. Parameters are stored for each patient in p#.mat, and for
every patient in par_all.mat.

Parameter estimation can be done for one_pop,two_pop,portz,and hirata models.

Created by Javier Baez Sep 2016. Modified by Tin Phan 2018.
%}
for yy = [2]
    % model = 'one_pop';                         % Choose which model to run fittings
    n = 3;                                       % number of periods to fit
    
    if yy == 1
        model = 'two_pop';
    elseif yy == 2
        model = 'newTwo_pop';
    end
    
    switch model                               % Number of parameters
        case 'two_pop'
            nParams = 20;
        case 'newTwo_pop'
            nParams = 17;
        otherwise
            warning('Unexpected model, choose a different one')
            return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%% Initialization of arrays and variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LIST = [1,6,14,15];%% done done
%     LIST = [17,19,28,29];%% done done
%     LIST = [32,36,37,39];%% done done
%     LIST = [51,52,54,55];%% done 
%     LIST = [60,62,63,66];%% done
%     LIST = [75,77,79,93];%% done
%     LIST = [105];%% done 
    LIST = [15];
      % Final patients: [1,6,14,15,17,19,25,28,29,32,36,37,39,51,52,54,55,60,62,63,66,75,77,79,93,100]; 
%%
      
%     LIST = [75 79 101];
%     LIST=[29];
    totalPatients = length(LIST);              % total number of patients
    patients = cell(1,totalPatients);          % creates a cell array to hold patient data
    counter = 1;                               % counter for patients and change array
    ind = cell(1,totalPatients);               % initializes the index for the S structure
    change = ones(1,n+1);                      % vector stores times when treatment changed
    par_store = zeros(nParams,totalPatients);  % used to store parameter values for each fit
%     options = optimset('Algorithm','interior-point','TolX',1e-16,'TolFun',1e-16,'TolCon'...
%     ,1e-13,'MaxIter',1000,'MaxFunEvals',1000);                     % Optimizer Options
    options = optimset('Algorithm','interior-point','TolX',1e-6,'TolFun',1e-6,'TolCon'...
    ,1e-6,'MaxIter',1000,'MaxFunEvals',1000);                     % Optimizer Options
    %%%%%%%%%%%%%%%%%%%% Creates Structure to store all patient data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = LIST
        patient = strcat('patient',num2str(ii));   % Patient with corresp number
        file = strcat('Data_62/',patient,'.txt');     % Complete name of file patient#.txt
        var = load(file);                          % holds variable just loaded
        patients(counter) = {var};                 % stores patient data into cell patients
        x = strcat('a',num2str(counter));          % Creates index for structure S
        ind{counter} = x;                          % puts all indexed in ind cell array
        counter = counter+ 1;                      % increases counter
    end
    S = cell2struct(patients,ind,2);              % cell to structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Runs the fitting for the selected patients %%%%%%%%%%%%%%%%%%%%%
    for i =1:totalPatients
        name = ['p',num2str(LIST(i))];    % where parameters will be saved
        index = char(ind(i));             % Index to calll specific patients
        patient = S.(index);              % calls specific patient from list S
        t = patient(:,2);                 % time vector in days
        psa = patient(:,3);               % psa vector of values
        androgen= patient(:,4);           % androgen vector of value

        treatment = patient(:,6);         % treatment vector of values
        %%%%%%%%%%%%%%%% Finds periods of on and off treatment %%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%% Bounds For Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch model
                case 'two_pop'
                a = max(androgen);
                b = min(androgen(change(1):change(4)));
                %  um            % q1             % q2              % c1
                LB(1) = 0.001;    LB(2) = b+.2;    LB(3) = 0.01;    LB(4) = 0.00001;
                UB(1) = 0.09;     UB(2) = b+.5;    UB(3) = b +.2;   UB(4) = 0.0001;
                %  K1            % b
                LB(5) = 0.8;      LB(6) = 0.0001;
                UB(5) = 1.7;      UB(6) = 0.1;
                %  sigma         % epsilon         % d1             % d2
                LB(7) = 0.001;   LB(8) = 0.0001;  LB(9) = 0.001;    LB(10) = 0.001;
                UB(7) = 1;       UB(8) = 1;       UB(9) = 0.09;     UB(10) = 0.01;
                %  R1            %  R2             % gamma1                 %  gamma2
                LB(11) = 0;      LB(12) = 1;       LB(13) = (0.08)*0.1;     LB(14) = 0.001;
                UB(11) = 3;      UB(12) = 6;       UB(13) = (0.08)*10;      UB(14) = 0.01;
                % dd1            % dd2             % Qm              % parameter u
                LB(15) = 1;      LB(16) = 1;       LB(17) = a - 5;   LB(20) = 0;
                UB(15) = 90;     UB(16) = 90;      UB(17) = a + 5;   UB(20) = 0;
                % X10             % X20                            
                LB(18) = 0.009;   LB(19) = 0.00001;       
                UB(18) = 0.02;    UB(19) = 0.001;    
                x0 = [androgen(1);0.01;0.0001;psa(1)]; %collects init cond in a column vector
                
            case 'newTwo_pop'
                a = max(androgen);
                b = min(androgen(change(1):change(4)));
                %  u1            % u2               % q1             % q2              
                LB(1) = 0.001;   LB(2) = 0.001;     LB(3) = b+.2;    LB(4) = 0.01;    
                UB(1) = 0.09;    UB(2) = 0.09;      UB(3) = b+.5;    UB(4) = b +.2;

                % b
                LB(5) = 0.0001;
                UB(5) = 0.1;
                % d1              
                LB(6) = 0.001;    
                UB(6) = 0.09;     
                %  R1                     
                LB(7) = 0;            
                UB(7) = 3;              
                % gamma1          
                LB(8) = (0.08)*0.1;      
                UB(8) = (0.08)*10;      
                % sigma          % epsilon
                LB(9) = 0.001;   LB(10) = 0.0001;
                UB(9) = 1;       UB(10) = 0.1;
                % X10             % X20                            
                LB(11) = 0.009;   LB(12) = 0.00001;       
                UB(11) = 0.02;    UB(12) = 0.001;       
                % A0                % u
                LB(13) = a-5;       LB(17) = 0;
                UB(13) = a+5;       UB(17) = 0;

                % d2              % R2           % m                          
                LB(14) = 0.01;    LB(15) = 1;    LB(16) = 0.01;     
                UB(14) = 0.09;    UB(15) = 6;    UB(16) = 0.9;   

                x0 = [androgen(1);0.01;0.0001;psa(1);0.9*androgen(1)]; %collects init cond in a column vector 
                %(A,x1,x2,P,Q)
        end
        %%%%%%%%%%%%%%%%%%%%%% Optimization Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IC = (LB+UB)/2; % initial parameter values
        [params,~] = fmincon(@(params)objective(params,psa,androgen,t,change,x0,n,model)...
            ,IC,[],[],[],[],LB,UB,[],options); % Finds optimized parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        [T1_test,Y1_test] = run_model(params,t,change,x0,n,model,androgen);
        figure(1);
        plot(t(1:change(n+1)),Y1_test(:,4)); hold on;
%         plot(T1_test,Y1_test(:,4)); hold on;
        scatter(t(1:change(n+1)),psa(1:change(n+1)))
        title(strcat('psa p',num2str(LIST(i))));
        savefig(strcat(int2str(patient(1,1)),model));
        hold off;
        figure(2);
        plot(t(1:change(n+1)),Y1_test(:,2),'r'); hold on;
        plot(t(1:change(n+1)),Y1_test(:,3),'b'); hold off;
        figure(3)
        plot(t(1:change(n+1)),Y1_test(:,1)); hold on;
        scatter(t(1:change(n+1)),androgen(1:change(n+1)))
        title(strcat('androgen A',num2str(LIST(i))));
        hold off;
           
        switch model
            case 'two_pop'
                x0 = [patient(4,1);params(18); params(19); patient(3,1)];
            case 'newTwo_pop'
                x0 = [patient(4,1);params(11); params(12); patient(3,1);0.9*patient(4,1)];
        end
        out(yy,i) = objective(params,psa,androgen,t,change,x0,n,model);
        
        par_store(:,i) = params; % stores parameters in a matrix to be used later for predictions and errors
        cd(strcat('parameters_',model))
        save(name,'params')      % saves parameters for individual patients in file p(patient#).mat
        cd ..
              
    end
    
    %%
    save(strcat('par_all_',model,'.mat'),'par_store') % saves the matrix of all patients parameter values
end
end

%%%%%%%%%%%%%%%%%  Minimizes the Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err]= objective(params,psadata,and_data,tdata,change,x0,n,model)
%{
 This function serves as the objective function that fmincon uses to find
 optimal parameter values.
 params  =  vector of parameters to fit
 psadata =  psa data
 change  =  vecotor that with the time steps at which treatmet is switched
 tdata   =  time data
 n       =  the number of periods of treatment to fit
%}
psadata = psadata(change(1):change(n+1));
and_data = and_data(change(1):change(n+1));

switch model

    case 'two_pop'
      x0 = [and_data(1);params(18);params(19);psadata(1)];
        [T1run,Y1run] = run_model(params,tdata,change,x0,n,model,and_data);        
        %
        train_first_cycle = tdata(change(n));
        timeF_first_cycle = tdata(1:change(n));

        train_second_cycle = tdata(change(n+1));
        timeF_second_cycle = tdata(change(n)+1:change(n+1));
        
        alpha_exponent_1 = 0.001;
        alpha_exponent_2 = 0.01;
        weight_first_cycle = 0.7;
        weight_second_cycle = 0.3;
       
        %
        PSA = Y1run(:,4);
        AND = Y1run(:,1);

        errp_first_cycle = sum((PSA(1:change(n))-psadata(1:change(n))).^2.*exp(-alpha_exponent_1*(ones(change(n),1)*train_first_cycle - timeF_first_cycle)))/length(PSA(1:change(n)));
        erra_first_cycle = sum((AND(1:change(n))-and_data(1:change(n))).^2.*exp(-alpha_exponent_1*(ones(change(n),1)*train_first_cycle - timeF_first_cycle)))/length(AND(1:change(n)));
        errp_second_cycle = sum((PSA(change(n)+1:change(n+1))-psadata(change(n)+1:change(n+1))).^2.*exp(-alpha_exponent_2*(ones(change(n+1)-change(n),1)*train_second_cycle - timeF_second_cycle)))/length(PSA(change(n)+1:change(n+1)));
        erra_second_cycle = sum((AND(change(n)+1:change(n+1))-and_data(change(n)+1:change(n+1))).^2.*exp(-alpha_exponent_2*(ones(change(n+1)-change(n),1)*train_second_cycle - timeF_second_cycle)))/length(AND(change(n)+1:change(n+1)));
        
        errp = weight_first_cycle*errp_first_cycle + weight_second_cycle*errp_second_cycle;
        erra = weight_first_cycle*erra_first_cycle + weight_second_cycle*erra_second_cycle;
        err = errp + erra;
        errorPSA = errp;
        errorAND = erra;
%         fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);        
            
    case 'newTwo_pop'
        x0 = [and_data(1);params(11);params(12);psadata(1);0.9*and_data(1)];
        [T1run,Y1run] = run_model(params,tdata,change,x0,n,model,and_data);        
        %
        train_first_cycle = tdata(change(n));
        timeF_first_cycle = tdata(1:change(n));

        train_second_cycle = tdata(change(n+1));
        timeF_second_cycle = tdata(change(n)+1:change(n+1));
        
        alpha_exponent_1 = 0.001;
        alpha_exponent_2 = 0.01;
        weight_first_cycle = 0.7;
        weight_second_cycle = 0.3;
       
        %
        PSA = Y1run(:,4);
        AND = Y1run(:,1);

        errp_first_cycle = sum((PSA(1:change(n))-psadata(1:change(n))).^2.*exp(-alpha_exponent_1*(ones(change(n),1)*train_first_cycle - timeF_first_cycle)))/length(PSA(1:change(n)));
        erra_first_cycle = sum((AND(1:change(n))-and_data(1:change(n))).^2.*exp(-alpha_exponent_1*(ones(change(n),1)*train_first_cycle - timeF_first_cycle)))/length(AND(1:change(n)));
        errp_second_cycle = sum((PSA(change(n)+1:change(n+1))-psadata(change(n)+1:change(n+1))).^2.*exp(-alpha_exponent_2*(ones(change(n+1)-change(n),1)*train_second_cycle - timeF_second_cycle)))/length(PSA(change(n)+1:change(n+1)));
        erra_second_cycle = sum((AND(change(n)+1:change(n+1))-and_data(change(n)+1:change(n+1))).^2.*exp(-alpha_exponent_2*(ones(change(n+1)-change(n),1)*train_second_cycle - timeF_second_cycle)))/length(AND(change(n)+1:change(n+1)));
        
        errp = weight_first_cycle*errp_first_cycle + weight_second_cycle*errp_second_cycle;
        erra = weight_first_cycle*erra_first_cycle + weight_second_cycle*erra_second_cycle;
        err = errp + erra;
        errorPSA = errp;
        errorAND = erra;
%
%         fprintf('PSA Error = %.4f \t Androgen Error = %.4f \n',errp, erra);
end
end
%%%%%%%%%%%%%%%%  Runs the Model and Generates synthetic data %%%%%%%%%%%%
function [T,y] = run_model(params,tdata,change,x0,n,model,androgen)
y = [];                                      % Initial vector for solution
T = [];
tint = tdata(change(n+1));                   % Time interval to run solution
for k = 1:n                                  % K number of cycles of treatment
    u = 1 - mod(k,2);
    switch model
        
        case 'two_pop'
            
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            
            if tint(end) >= tdata(change(k))
                [trun,Yrun]=ode15s(@(t,x) two_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0,options);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4)];
                if k < n
                    y = [y; Yrun(1:end-1,:)];
                    T = [T; trun(1:end-1,:)];
                elseif k == n
                    y = [y; Yrun(1:end,:)]; 
                    T = [T; trun(1:end-1,:)];
                end
        end
        
        %%%%%%%%%%%
        case 'newTwo_pop'
            
            options = odeset('RelTol',1e-12,'AbsTol',1e-14);
            
            if tint(end) >= tdata(change(k))
              
                [trun,Yrun]=ode15s(@(t,x) newTwo_pop(t,x,[params(1:end-1),u]),tdata(change(k):change(k+1)),x0,options);%,options);
                x0 = [Yrun(end,1);Yrun(end,2);Yrun(end,3);Yrun(end,4);Yrun(end,5)];
                if k < n
                    y = [y; Yrun(1:end-1,:)]; %#ok<AGROW>
                    T = [T; trun(1:end-1,:)]; %#ok<AGROW>
                elseif k == n
                    y = [y; Yrun(1:end,:)];   %#ok<AGROW>
                    T = [T; trun(1:end,:)]; %#ok<AGROW>
                end
            end
    end
end
end
% model ODE functions

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
