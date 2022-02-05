%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNessterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Nonstrongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

clear all

% global variables
global delta M gamma lambda c_0 c_10 r tauMin tauMax tauMed c zeta cTilde_0 cTilde_10 d_0 d_10 alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
M = 2;
zeta = 2;

% Heavy Ball constants
gamma = 2/3; % 
lambda = 200; % 40

% Uniting parameters for \mathcal{U}_0 and \mathcal{T}_{1,0}:
c_0 = 2000;  
c_10 = 1154.148149;

% These are the same, since L = x^2
alpha = 1;

% eps_0 has to be bigger than eps_10
eps_0 = 10;
eps_10 = 5; % 15

cTilde_0 = eps_0*alpha
cTilde_10 = eps_10*alpha
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)

tauMin = (1+sqrt(7))/2;
c = 0.5;
deltaMed = 8112; 
r = 21; 

delta = 0.2;

% initial conditions:
z1_0 = 20; 
z2_0 = 0;
z2_00 = 20;
q_0 = 1;
tau_0 = 0;
tauPN_0 = tauMin; 

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN_HBF=[0 2000];
TSPAN=[0 100];
JSPAN = [0 20000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

%% Simulating Uniting and HAND-1 from z1(0,0) = 20:

% Uniting
[tU9,jU9,xU9] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU9 = CalcLVector(xU9,tU9);
deltaVecU9 = timeToConv(xU9,tU9);

% HAND-1
[tH9,jH9,xH9] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH9 = CalcLVector(xH9,tH9);
deltaVecH9 = timeToConv(xH9,tH9);

% heavy ball
[tHBF9,jHBF9,xHBF9] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF9 = timeToConv(xHBF9,tHBF9);

% Nesterov
[tNest9,jNest9,xNest9] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest9 = timeToConv(xNest9,tNest9);

%% Simulating Uniting and HAND-1 from z1(0,0) = 30:
% Retune r and deltaMed for HAND-1
deltaMed = 18110;
r = 31;

c_0 = 3000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 2503.083335; % Too high: 2503.08334  Too low: 2503.08333
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.3;

% initial conditions
z1_0 = 30; 
z2_00 = 30;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU8,jU8,xU8] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU8 = CalcLVector(xU8,tU8);
deltaVecU8 = timeToConv(xU8,tU8);

% HAND-1
[tH8,jH8,xH8] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH8 = CalcLVector(xH8,tH8);
deltaVecH8 = timeToConv(xH8,tH8);

% heavy ball
[tHBF8,jHBF8,xHBF8] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF8 = timeToConv(xHBF8,tHBF8);

% Nesterov
[tNest8,jNest8,xNest8] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest8 = timeToConv(xNest8,tNest8);

%% Simulating Uniting and HAND-1 from z1(0,0) = 40:
% Retune deltaMed for HAND-1
deltaMed = 32075; 
r = 41;

c_0 = 5000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 4391.5925953; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.4;

% initial conditions
z1_0 = 40; 
z2_00 = 40;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU7,jU7,xU7] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU7 = CalcLVector(xU7,tU7);
deltaVecU7 = timeToConv(xU7,tU7);

% HAND-1
[tH7,jH7,xH7] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH7 = CalcLVector(xH7,tH7);
deltaVecH7 = timeToConv(xH7,tH7);

% heavy ball
[tHBF7,jHBF7,xHBF7] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF7 = timeToConv(xHBF7,tHBF7);

% Nesterov
[tNest7,jNest7,xNest7] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest7 = timeToConv(xNest7,tNest7);

%% Simulating Uniting and HAND-1 from z1(0,0) = 50:
% Retune deltaMed for HAND-1
deltaMed = 50000;

r = 51;

c_0 = 7000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta for Uniting:
c_10 = 6819.67593; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.5;

% initial conditions
z1_0 = 50; 
z2_00 = 50;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU6,jU6,xU6] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU6 = CalcLVector(xU6,tU6);
deltaVecU6 = timeToConv(xU6,tU6);

% HAND-1
[tH6,jH6,xH6] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH6 = CalcLVector(xH6,tH6);
deltaVecH6 = timeToConv(xH6,tH6);

% heavy ball
[tHBF6,jHBF6,xHBF6] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF6 = timeToConv(xHBF6,tHBF6);

% Nesterov
[tNest6,jNest6,xNest6] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest6 = timeToConv(xNest6,tNest6);

%% Simulating Uniting and HAND-1 from z1(0,0) = 60:
% Retune deltaMed for HAND-1
deltaMed = 71875; 
r = 61;

c_0 = 10500;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 9787.33334; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.6;

% initial conditions
z1_0 = 60; 
z2_00 = 60;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU5,jU5,xU5] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU5 = CalcLVector(xU5,tU5);
deltaVecU5 = timeToConv(xU5,tU5);

% HAND-1
[tH5,jH5,xH5] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH5 = CalcLVector(xH5,tH5);
deltaVecH5 = timeToConv(xH5,tH5);

% heavy ball
[tHBF5,jHBF5,xHBF5] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF5 = timeToConv(xHBF5,tHBF5);

% Nesterov
[tNest5,jNest5,xNest5] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest5 = timeToConv(xNest5,tNest5);

%% Simulating Uniting and HAND-1 from z1(0,0) = 70:
% Retune deltaMed for HAND-1
deltaMed = 97700; 
r = 71;

c_0 = 14000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 13294.564825; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.7;

% initial conditions
z1_0 = 70; 
z2_00 = 70;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU4,jU4,xU4] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU4 = CalcLVector(xU4,tU4);
deltaVecU4 = timeToConv(xU4,tU4);

% HAND-1
[tH4,jH4,xH4] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH4 = CalcLVector(xH4,tH4);
deltaVecH4 = timeToConv(xH4,tH4);

% heavy ball
[tHBF4,jHBF4,xHBF4] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF4 = timeToConv(xHBF4,tHBF4);

% Nesterov
[tNest4,jNest4,xNest4] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest4 = timeToConv(xNest4,tNest4);

%% Simulating Uniting and HAND-1 from z1(0,0) = 80:
% Retune deltaMed for HAND-1
deltaMed = 127550; 
r = 81;

c_0 = 18000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 17341.370385; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.8;

% initial conditions
z1_0 = 80; 
z2_00 = 80;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU3,jU3,xU3] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU3 = CalcLVector(xU3,tU3);
deltaVecU3 = timeToConv(xU3,tU3);

% HAND-1
[tH3,jH3,xH3] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH3 = CalcLVector(xH3,tH3);
deltaVecH3 = timeToConv(xH3,tH3);

% heavy ball
[tHBF3,jHBF3,xHBF3] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF3 = timeToConv(xHBF3,tHBF3);

% Nesterov
[tNest3,jNest3,xNest3] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest3 = timeToConv(xNest3,tNest3);

%% Simulating Uniting and HAND-1 from z1(0,0) = 90:
% Retune deltaMed for HAND-1
deltaMed = 161300; 
r = 91;

c_0 = 23000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 21927.75002;  
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.9;

% initial conditions
z1_0 = 90; 
z2_00 = 90;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU2,jU2,xU2] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU2 = CalcLVector(xU2,tU2);
deltaVecU2 = timeToConv(xU2,tU2);

% HAND-1
[tH2,jH2,xH2] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH2 = CalcLVector(xH2,tH2);
deltaVecH2 = timeToConv(xH2,tH2);

% heavy ball
[tHBF2,jHBF2,xHBF2] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF2 = timeToConv(xHBF2,tHBF2);

% Nesterov
[tNest2,jNest2,xNest2] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest2 = timeToConv(xNest2,tNest2);

%% Simulating Uniting and HAND-1 from z1(0,0) = 100:
% Retune deltaMed for HAND-1
deltaMed = 199000; 
r = 101;

c_0 = 28000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 27053.703733; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 1;

% initial conditions
z1_0 = 100; 
z2_00 = 100;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU1,jU1,xU1] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU1 = CalcLVector(xU1,tU1);
deltaVecU1 = timeToConv(xU1,tU1);

% HAND-1
[tH1,jH1,xH1] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH1 = CalcLVector(xH1,tH1);
deltaVecH1 = timeToConv(xH1,tH1);

% heavy ball
[tHBF1,jHBF1,xHBF1] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF1 = timeToConv(xHBF1,tHBF1);

% Nesterov
[tNest1,jNest1,xNest1] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest1 = timeToConv(xNest1,tNest1);

%% Simulating Uniting and HAND-1 from z1(0,0) = 110:
% Retune deltaMed for HAND-1
deltaMed = 240700; 
r = 111;

c_0 = 34000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 32719.23153;  % Too low: 32719.2315 Too high: 32719.23154
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 1.1;

% initial conditions
z1_0 = 110; 
z2_00 = 110;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting
[tU,jU,xU] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
lU = CalcLVector(xU,tU);
deltaVecU = timeToConv(xU,tU);

% HAND-1
[tH,jH,xH] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');
lH = CalcLVector(xH,tH);
deltaVecH = timeToConv(xH,tH);

% heavy ball
[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');
deltaVecHBF = timeToConv(xHBF,tHBF);

% Nesterov
[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');
deltaVecNest = timeToConv(xNest,tNest);

%% Group together the trajectories, to plot:
minarc = min([length(xU),length(xU1),length(xU2),length(xU3),length(xU4),length(xU5),length(xU6),length(xU7),length(xU8),length(xU9)]); 
ta = [tU(1:minarc),tU1(1:minarc),tU2(1:minarc),tU3(1:minarc),tU4(1:minarc),tU5(1:minarc),tU6(1:minarc),tU7(1:minarc),tU8(1:minarc),tU9(1:minarc)]; 
ja = [jU(1:minarc),jU1(1:minarc),jU2(1:minarc),jU3(1:minarc),jU4(1:minarc),jU5(1:minarc),jU6(1:minarc),jU7(1:minarc),jU8(1:minarc),jU9(1:minarc)]; 
xa = [xU(1:minarc,1),xU1(1:minarc,1),xU2(1:minarc,1),xU3(1:minarc,1),xU4(1:minarc,1),xU5(1:minarc,1),xU6(1:minarc,1),xU7(1:minarc,1),xU8(1:minarc,1),xU9(1:minarc,1)]; 
la = [lU(1:minarc).',lU1(1:minarc).',lU2(1:minarc).',lU3(1:minarc).',lU4(1:minarc).',lU5(1:minarc).',lU6(1:minarc).',lU7(1:minarc).',lU8(1:minarc).',lU9(1:minarc).'];

minarc = min([length(xH),length(xH1),length(xH2),length(xH3),length(xH4),length(xH5),length(xH6),length(xH7),length(xH8),length(xH9)]); 
tb = [tH(1:minarc),tH1(1:minarc),tH2(1:minarc),tH3(1:minarc),tH4(1:minarc),tH5(1:minarc),tH6(1:minarc),tH7(1:minarc),tH8(1:minarc),tH9(1:minarc)];
jb = [jH(1:minarc),jH1(1:minarc),jH2(1:minarc),jH3(1:minarc),jH4(1:minarc),jH5(1:minarc),jH6(1:minarc),jH7(1:minarc),jH8(1:minarc),jH9(1:minarc)];
xb = [xH(1:minarc,1),xH1(1:minarc,1),xH2(1:minarc,1),xH3(1:minarc,1),xH4(1:minarc,1),xH5(1:minarc,1),xH6(1:minarc,1),xH7(1:minarc,1),xH8(1:minarc,1),xH9(1:minarc,1)];
lb = [lH(1:minarc).',lH1(1:minarc).',lH2(1:minarc).',lH3(1:minarc).',lH4(1:minarc).',lH5(1:minarc).',lH6(1:minarc).',lH7(1:minarc).',lH8(1:minarc).',lH9(1:minarc).'];

%% Plotting
figure(1) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
subplot(1,2,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
axis([0 3 -150 150])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
subplot(1,2,2), plotHarc(tb,jb,xb,[],modificatorF,modificatorJ);
axis([0 3 -150 150])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
saveas(gcf,'Plots\TrajectoryPlotsNSC','png')

figure(2)
clf
subplot(1,2,1), semilogy(ta,la,'LineWidth',1.5);
axis([0 3 10^(-20) 10^5])
legend({'z_1(0,0)=110','z_1(0,0)=100','z_1(0,0)=90','z_1(0,0)=80','z_1(0,0)=70','z_1(0,0)=60','z_1(0,0)=50','z_1(0,0)=40','z_1(0,0)=30','z_1(0,0)=20'},'Location','northeast')
grid on
ylabel('L(z_1)-L^*','Fontsize',16)
xlabel('t','Fontsize',16)
subplot(1,2,2), semilogy(tb,lb,'LineWidth',1.5);
axis([0 3 10^(-20) 10^5])
legend({'z_1(0,0)=110','z_1(0,0)=100','z_1(0,0)=90','z_1(0,0)=80','z_1(0,0)=70','z_1(0,0)=60','z_1(0,0)=50','z_1(0,0)=40','z_1(0,0)=30','z_1(0,0)=20'},'Location','southwest')
grid on
ylabel('L(z_1)-L^*','Fontsize',16)
xlabel('t','Fontsize',16)
saveas(gcf,'Plots\TrajectoryPlotsNSCLValue','png')