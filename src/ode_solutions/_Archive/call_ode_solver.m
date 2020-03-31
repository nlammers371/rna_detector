% Numerical simulations to simple ODEs to convey general pricniples behind
% NCR
clear
close all
addpath('../utilities')


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

x0_vec = [A0 S0 AS0 P0];

% rate parameters 
param_struct = struct;

param_struct.koff1 = 9;%1e-3;
param_struct.kon1 = 1;%1e-3;
param_struct.kc1 = 1;

% try solver
tic
[t15,y15] = ode15s(@ncr_damper_v1,[0 360],x0_vec);
toc

% [t,y] = ode45(@odefun_ncr_v1,[0 total_time],0)
