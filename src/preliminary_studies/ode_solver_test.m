% Numerical simulations to simple ODEs to convey general pricniples behind
% NCR
clear
close all

addpath('utilities')
FigPath = '../fig/preliminary_studies/';
mkdir(FigPath)

DataPath = '../out/preliminary_studies/';
mkdir(DataPath)


% NCR
total_time = 3600;

% initial concentration of active complex 
% A0_vec = [1e-8 1e-6 1e-4 1e-2 1e0]*10^-9; 

% define valency parameter
valency_vec = [0, 1];

% other initial concentrations (in molar)
S0 = 2e4; % dark reporter
P0 = 0; % cleaved (fluorescent) reporter
B0 = 2e4; % caged amplifier
AB0 = 0; % cage-active complex
AS0 = 0;
A0 = 100;

x0_vec = [A0 S0 AS0 P0 B0 AB0];


dt = 0.01;%1/max([kon1 kon2]);
n_steps = round(total_time/dt);


% try solver
tic
[t15,y15] = ode15s(@odefun_ncr_v1,[0 100],x0_vec);
toc
tic
[t45,y45] = ode45(@odefun_ncr_v1,[0 100],x0_vec);
toc
% [t,y] = ode45(@odefun_ncr_v1,[0 total_time],0)
