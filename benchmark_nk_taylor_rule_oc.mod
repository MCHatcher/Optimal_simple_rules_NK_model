//Benchmark 3-equation New Keynesian model: Taylor rules vs optimal commitment

//The basic idea is to simulate both cases using the same set of shocks
//Hence, two different model economies are built within one mod file
//The optimal Taylor rule coefficients on inflation and the output gap are computed using the loops facility in Dynare 
 
// Parameter values are taken from Chapter 5 in Gali(2008)
// Written by Michael Hatcher (Southampton), building on the benchmark code of Ding Liu (SWUFE)

var x //welfare-relevant output gap
    pi //inflation
    r //nominal interest rate
    x_oc //gap under optimal commitment policy
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
r = phi_pi*pi + phi_x*x;
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

init_coef_pi = 1.01;
ncoefs_pi = 61; //number of coefficients in inflation direction
ncoefs_x = 61; //number of coefficients in output gap direction
max_pi_coef = 15; 
max_x_coef = 10;

welfare_loss = zeros(ncoefs_pi,ncoefs_x);
welfare_loss_oc = zeros(ncoefs_pi,ncoefs_x);

for j=1:ncoefs_pi

    for k=1:ncoefs_x

            coef(j) = init_coef_pi + (max_pi_coef-init_coef_pi)*(j-1)/ncoefs_pi;
            phi_pi = coef(j);

            coef1(k) = (k-1)*max_x_coef/ncoefs_x;
            phi_x = coef1(k);

options_.qz_criterium = 1+1e-6;
steady;
check;
stoch_simul(order=1, periods=0, irf=0, noprint); //periods=0: theoretical moments option 
//stoch_simul(order=1, periods=11100, drop=100, irf=0, noprint); //simulated moments option (takes several hours) 
var_x(j,k) = oo_.var(1,1); % output gap variance
var_pi(j,k) = oo_.var(2,2); % inflation variance
welfare_loss(j,k) = -(var_pi(j,k)+lambda_y*var_x(j,k));

var_x_oc(j,k) = oo_.var(4,4); % output gap variance
var_pi_oc(j,k) = oo_.var(5,5); % inflation variance
welfare_loss_oc(j,k) = -(var_pi_oc(j,k)+lambda_y*var_x_oc(j,k));

end;
end;

welfare_loss;

//Optimal output gap coefficient
MN = min(abs(welfare_loss));  //row vector containing min for each column
[MN1, Index_x] = min(MN);  //finds which row has lowest loss and records location

//Optimal inflation coefficient
MN2 = min(abs(welfare_loss'));  //row vector containing min for each column (note: transposed)
[MN3, Index_pi] = min(MN2); //finds which row has lowest loss and records location (note: tranposed)

Index_x;     //Index for optimal coefficient on output gap
Index_pi;    //index for optimal coefficient on inflation

Optimal_pi_coef = init_coef_pi + (max_pi_coef-init_coef_pi)*(Index_pi-1)/ncoefs_pi  //Loss-minmising inflation coefficient in Taylor rule
Optimal_x_coef = (Index_x-1)*max_x_coef/ncoefs_x  //Loss-minmising gap coefficient in Taylor rule

Min_loss_rule = welfare_loss(Index_pi,Index_x)
Loss_oc = welfare_loss_oc(Index_pi,Index_x)

figure(1)
subplot(1,2,1), surf(coef1, coef, welfare_loss), title('Social loss under Taylor rule'), 
ylabel('pi coef'), xlabel('x coef')
subplot(1,2,2), surf(coef1, coef, welfare_loss_oc), 
ylabel('pi coef'), xlabel('x coef');
title('Social loss under optimal commitment')

figure(2)
subplot(1,2,1), surf(coef1, coef, var_x), title('Output gap variance under Taylor rule'), 
ylabel('pi coef'), xlabel('x coef')
subplot(1,2,2),surf(coef1, coef, var_x_oc), title('Output gap variance under optimal commitment'), 
ylabel('pi coef'), xlabel('x coef')

figure(3)
subplot(1,2,1), surf(coef1, coef, var_pi), title('Inflation  variance under Taylor rule'), 
ylabel('pi coef'), xlabel('x coef')
subplot(1,2,2),surf(coef1, coef, var_pi_oc), title('Inflation variance under optimal commitment'), 
ylabel('pi coef'), xlabel('x coef')
