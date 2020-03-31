% Script to explore behavior of simplified model that uses MM kinetics to
% for all first order interactions

clear
close all
addpath('../utilities')

FigPath = '../fig/ode_studies/';
mkdir(FigPath)
DataPath = '../out/ode_studies/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key (first order)
%   species 1: A13I (caged activator)
%   species 2: A13 (free activator and/or virus)
%   species 3: C13 (free Cas13)
%   species 4: S (dark reporter)
%   species 5: F (Cleaved (fluorescent) reporter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first define a set of ODEs for association/disassociation that we will
% assume to be in rapid equilibrium wrpt cleavage activity

% reactants (first order)
reactant_names = {'A13' 'C13' 'S' 'A13I' 'A13_C13' 'A13I_C13' 'S_C13' 'S_A13_C13' 'A13I_A13_C13'};
ractant_vec = sym(reactant_names); 


% rates
rate_names = {'kon', 'koff_s', 'koff_ns', 'kcat_low', 'kcat_high'}; 
% dydt_sym = {};
dydt_fun = {};

% A
eqA = -kon*ractant_vec(1)*ractant_vec(2) + koff_s*ractant_vec(5);
dydt_fun{1} = @(y) -kon*y(1)*y(2) + koff_s*y(5);
% C
eqC = -kon*C13*(A13 + S + A13I) + koff_s*A13_C13 + koff_ns*(S+A13I) == 0;
% S
eqS = -kon*S*(C13 + A13_C13) + koff_ns*(S_C13 + S_A13_C13) == 0;
% A13I
eqAI = -kon*A13I*(C13 + A13_C13) + koff_ns*(A13I_C13 + A13I_A13_C13) == 0;
% CA
eqCA = kon*C13*A13 - koff_s*A13_C13 -kon*A13_C13*(S + A13I) == 0;
% C:AI
eqCAI = kon*C13*A13I - koff_ns*A13I_C13 == 0;
% C:S
eqCS = kon*C13*S - koff_ns*S_C13 == 0;
% CA:AI
eqCAAI = kon*A13_C13*A13I - koff_ns*A13I_A13_C13 == 0;
% CA:S
eqCAS = kon*A13_C13*S - koff_ns*S_A13_C13 == 0;
% 
% % define system of equations and corresponding variables
% eq_sys = [eqA eqC eqS eqAI eqCA eqCAI eqCS eqCAAI eqCAS];
% sys_vars = [A13 C13 S A13I  A13_C13 A13I_C13 S_C13 A13I_A13_C13 S_A13_C13];
% % solve system
% sys_sol = solve(eqA,A13)







% define system of equations and corresponding variables
eq_sys = [eqCA eqC eqA eqS eqCS];
sys_vars = [A13_C13 C13 A13 S S_C13];

% solve system
% x0 = rand(1,5);
tic
sys_sol = solve(eq_sys,sys_vars);
toc

% check that results make sense. Lets titrate A while holding S and C fixed

