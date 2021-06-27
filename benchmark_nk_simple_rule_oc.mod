//Benchmark 3-equation New Keynesian model: simple rule vs optimal commitment

//The basic idea is to simulate both cases using the same set of shocks
//Hence, two different model economies are built within one mod file
//The optimal feedback coefficient on inflation in the simple rule is computed using the loops facility in Dynare 
 
// Parameter values are taken from Chapter 5 in Gali(2008)
// Written by Michael Hatcher (Southampton), building on the benchmark code of Ding Liu (SWUFE)

var x //welfare-relevant output gap
    pi //inflation
    r //nominal interest rate
    x_oc //output gap under optimal commitment policy
    pi_oc //inflation under optimal commitment policy
    r_oc //nominal interest rate under optimal commitment policy
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
rho_u = 0.5;
rho = 0.5;

model(linear);
//1. IS equation
x = x(+1)-sigma*(r-pi(+1)-r_n);
//2. NK Phillips curve
pi = beta*pi(+1)+kappa*x+u;
//3. Interest rate rule
r = phi_pi*pi;
//4. cost-push shock
u = rho_u*u(-1)+e_u;
//5. IS equation (optimal commitment)
x_oc = x_oc(+1)-sigma*(r_oc-pi_oc(+1)-r_n);
//6. NK Phillips curve (optimal commitment) 
pi_oc = beta*pi_oc(+1)+kappa*x_oc+u;
//7. Optimal commitment policy
pi_oc = -lambda_y/kappa*(x_oc-x_oc(-1));
//8. Natural rate shock
r_n = rho*r_n(-1) + e_rn;
end;

steady_state_model;
x=0;
r=0;
pi=0;
x_oc=0;
r_oc=0;
pi_oc=0;
u=0;
r_n=0;
end;

shocks;
var e_rn;  stderr 0.5;
var e_u;  stderr 1;
end;

close all;

init_coef = 1.3;
ncoefs = 1100; //number of inflation coefficients in loop
max_coef = 12;

for j=1:ncoefs

coef(j) = init_coef + (max_coef-init_coef)*(j-1)/ncoefs;
phi_pi = coef(j);

options_.qz_criterium = 1+1e-6;
steady;
check;
stoch_simul(order=1, periods=0, irf=0, noprint); //periods=0: theoretical moments option 
//stoch_simul(order=1, periods=11100, drop=100, irf=0, noprint); //simulated moments option (takes several minutes) 
var_x(j) = oo_.var(1,1); % output gap variance
var_pi(j) = oo_.var(2,2); % inflation variance
welfare_loss(j) = -(var_pi(j)+lambda_y*var_x(j));

var_x_oc(j) = oo_.var(4,4); % output gap variance
var_pi_oc(j) = oo_.var(5,5); % inflation variance
welfare_loss_oc(j) = -(var_pi_oc(j)+lambda_y*var_x_oc(j));

end;

[Min_Loss_Rule, Index] = min(abs(welfare_loss)); //Minimum value and location of min social loss (simple rule)
Optimal_pi_coef = init_coef + (max_coef-init_coef)*(Index-1)/ncoefs  //Loss-minmising inflation coefficient in simple rule

Min_Loss_Rule = -Min_Loss_Rule

figure(1)
hold on,
plot(coef, welfare_loss,'b'),
plot(coef, welfare_loss_oc,'r'), xlabel('Inflation reaction coefficient'), ylabel('Welfare loss');
title('Social loss: rule (blue) vs optimal commitment (red)')

figure(2)
hold on,
plot(coef, var_x,'b'), 
plot(coef, var_x_oc,'r'), xlabel('Inflation reaction coefficient'), ylabel('Output gap variance');
title('Output gap variance: rule (blue) vs optimal commitment (red)')

figure(3)
hold on,
plot(coef, var_pi,'b'), 
plot(coef, var_pi_oc,'r'), xlabel('Inflation reaction coefficient'), ylabel('Inflation variance');
title('Inflation variance: rule (blue) vs optimal commitment (red)')


