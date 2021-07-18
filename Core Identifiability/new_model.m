% model ODE function
function dxdt = new_model(t,x,p) 
mu = p(1);  q2= p(2); gamma1 = p(4); d1 = p(3); A0 = p(5); 

dd1 = 5; dd2 = 5;
q1= 0.5;
b = 0.05; R1 = 3; 
sigma = 0.5; epsilon = 0.05;
d2 = 0.009; R2 = 4; m = 0.5;
c = 0.00015; K = 1; gamma2 = 0.005;
%separates solutions 
Q = x(5); X = x(2); Y = x(3); P = x(4); A = x(1);
%%%%%%%%%%%%%%%%%%%  Parameters for ODE System %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q1>Q
    ux = 0;
else
    ux = mu*(1 - abs(q1)/Q);
end
if q2 > Q
    uy = 0;
else
    uy = mu*(1 - abs(q2)/Q);
end
Dx = d1*R1/(Q+R1);    Dy = d2*R2/(Q+R2);
mxy = c*K/(Q + K);
%%%%%%%%%%%%%%%%%%% ODE system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dX = ux*X - Dx*X - dd1*X*X - mxy*X;
dY = uy*Y - Dy*Y - dd2*Y*Y + mxy*X;
dA = gamma1*(A0 - A) - gamma1*A0*(1-u(t)) + gamma2;
dQ = m*(A - Q) - (ux*Q*X + uy*Q*Y)/(X+Y);
dP = b*Q + sigma*(Y*Q + X*Q) - epsilon*P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxdt = [dA; dX; dY; dP; dQ]; %Puts the ode in a column vector
end
