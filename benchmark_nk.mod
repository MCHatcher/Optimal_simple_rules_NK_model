//Benchmark 3-equation New Keynesian model with Taylor rule
 
// Parameter values are taken from Chapter 5 in Gali(2008)
// Written by Ding Liu (SWUFE) and Michael Hatcher (Southampton)

var x //welfare-relevant output gap
    pi //inflation
    r //nominal interest rate
    r_n //natural rate shock
    u; // cost-push shock in (3) p.97

varexo e_rn //innovation of natural rate shock
       e_u; // innovation of cost-push shock

parameters 
    beta //discount factor 
    alpha //capital share
    varphi //Frisch elasticity
    theta //Calvo parameter
    sigma //Risk aversion
    epsilon //Elasticity of substitution
    phi_pi //Taylor rule feedback inflation
    phi_x //Taylor rule feedback output gap
    rho_u //Autocorrelation of cost-push shock
    rho //Autocorrelation of natural rate shock
    lambda_y //Weight of output gap in the loss function
    lambda 
    kappa; //Slope of NK Phillips curve

beta = 0.99;
alpha = 0;
theta = 0.75;
epsilon = 10;
sigma = 1; 
varphi = 1;
lambda  = (1-theta)*(1-beta*theta)/theta*(1-alpha)/(1-alpha+alpha*epsilon); // p.47
kappa  = lambda*(sigma+(varphi+alpha)/(1-alpha)); // p.49
lambda_y = kappa/epsilon; //p.96
phi_pi = 1.5; // rule coef on inflation
phi_x = 0.5; // rule coef on welfare-relevant output gap
rho_u = 0.5;
rho = 0.5;

model(linear);
//1. IS equation
x = x(+1)-sigma*(r-pi(+1)-r_n);
//2. NK Phillips curve
pi = beta*pi(+1)+kappa*x+u;
//3. Interest rate rule
r = phi_pi*pi+phi_x*x;
//4. cost-push shock
u = rho_u*u(-1)+e_u;
//5. Natural rate shock
r_n = rho*r_n(-1) + e_rn;
end;

steady_state_model;
x=0;
r=0;
pi=0;
u=0;
r_n=0;
end;

shocks;
var e_rn;  stderr 1;
var e_u;  stderr 1;
end;

options_.qz_criterium = 1+1e-6;

steady;
check;
stoch_simul(periods = 1000, irf = 20);
