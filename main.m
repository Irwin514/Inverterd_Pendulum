clear all 
close all 
clc
pkg load control



% Parameters Initialization
m1 = .57;           % mass of the cart
m2 = 0.23;          % mass of the Rod
l = 0.6413;         % length of the Rod
g = 9.81;           % gravity Earth
d1 = 1;             % Damping Term on the cart 

% Jacobians Matrix:
A = [0 1 0 0;
    0 -d1/m1  m2*g/m1 0;
    0 0 0 1;
    0 -d1/(m1*l) (m1+m2)*g/(m1*l) 0];
    
B = [0; 1/m1; 0; 1/(m1*l)];

%Compute eigenvalues
eig(A);

% check rank to see if the system is controllable
rank(ctrb(A,B));

% Q matrix for Optimization Problem
Q = [10 0 0 0;   % x_ = [x x_d theta theta_d]
      0 100 0 0;
      0 0 10000 0;
      0 0 0 1];
R = 0.00001;

% Solving Algebraic Riccati Equation ODE to determine K (Gain)
K = lqr(A,B,Q,R);


% create time vector
tspan = 0:0.1:15;

% initial condition: [x, x_d, theta, theta_d]
x0 = [0; 0; -pi; 2]; 
%x0 = [-1; 0; 0; 0.5];

% ode options
odeset('relTol', 1e-7, 'absTol', realmin);

% Solving the Dormand_Prince Method of Order 4
[t_out,x_out] = ode45( @(t,x)evaluateOde(t,x,m1,m2,l,g,d1,K), tspan, x0, odeset);

% draw the pendulum animation
drawPendulum(t_out, x_out, m1, m2, g, l);