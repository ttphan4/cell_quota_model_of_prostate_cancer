function err = Objective(par,psa_data,and_data,tdata,x0)
opts = odeset('RelTol',1e-6,'AbsTol',1e-10);
[t,est] = ode23s(@new_model,tdata,x0,opts,par);
err = sum((est(:,1)-and_data).^2+(est(:,4)-psa_data).^2);
e = err;
end