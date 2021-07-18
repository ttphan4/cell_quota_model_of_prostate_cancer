function err = Objective_fixed_mu(par_fixed_mu,psa_data,and_data,tdata,x0)
global params e fixed_mu
params = [fixed_mu,par_fixed_mu];
opts = odeset('RelTol',1e-6,'AbsTol',1e-10);
[t,est] = ode23s(@new_model_fixed_mu,tdata,x0,opts);
err = sum((est(:,1)-and_data).^2+(est(:,4)-psa_data).^2);
e = err;
end