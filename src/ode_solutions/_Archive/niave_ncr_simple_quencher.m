% Script to illustrate the sensitivity issue with naive NCR

% Even "negative" samples with low level background activity will trigger
% the cascade, and the timing with which the signals for positve and
% negative samples is quite close, even given relatively large initial
% differences

clear
close all
addpath('../utilities')

FigPath = '../fig/ode_studies/';
mkdir(FigPath)
DataPath = '../out/ode_studies/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key:   
%   species 1: AI (caged activator)
%   species 2: A (free activator)
%   species 3: CA (Cas13:Activator)    
%   species 4: CAAI (Cas13:Activator:Caged activator)
%   species 5: C (free Cas13)
%   species 6: S (dark reporter)
%   species 7: SCA (dark reporter  + Cas13:activator)
%   species 8: F (Cleaved (fluorescent) reporter)

%   species 9: QI ( caged quencher) 
%   species 10: Q (free quencher)
%   species 11: QA (quenched activator)
%   species 12: CAQI (Cas13:Activator:Caged Quencher)

%   species 13: N (background nuclease contaminant)
%   species 14: NAI (bkg nuclease : Caged activator)
%   species 15: NQI (bkg nuclease : Caged quencher)
%   species 16: NS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reaction key:    
%    reaction 1: C + A -> CA    %kon
%    reaction 2: CA -> C + A    %koff

%    reaction 3: CA + AI -> CAAI    %kon   

%    reaction 4: CAAI -> CA + AI    %koff        
%    reaction 5: CAAI -> CA + A    %kc

%    reaction 6: S + CA -> SCA    %kon
%    reaction 7: SCA -> S + CA      %koff 
%    reaction 8: SCA -> CA + F     %kcat   

%    reaction 9: CA + QI -> CAQI    %kon
%    reaction 10: CAQI ->  CA + QI   %koff
%    reaction 11: CAQI ->  CA + Q   %kcat_ns

%    reaction 12: Q + A ->  QA   %kon
%    reaction 13: QA -> Q + A   %koff

%    reaction 14: N + AI -> NAI    %kon   
%    reaction 15: NAI -> N + AI    %koff        
%    reaction 16: NAI -> N + A    %kcat

%    reaction 17: N + QI -> NQI    %kon   
%    reaction 18: NQI -> N + QI    %koff        
%    reaction 19: NQI -> N + Q    %kcat

%    reaction 20: S + N -> NS %kon
%    reaction 21: SN -> N + S %koff
%    reaction 22: SN -> N + F %kcat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define stoichiometry matrix
Q = zeros(16,22);    

Q(1,[3 4 14 15]) = [-1 1 -1 1]; % AI
Q(2,[1 2 5 12:13]) = [-1 1 1 -1 1]; %A
Q(3,1:11) = [1 -1 -1 1 1 -1 1 1 -1 1 1]; % CA
Q(4,[3 4 5]) = [1 -1 -1]; %CAAI
Q(5,[1 2]) = [-1 1]; %C
Q(6,[6 7 20 21]) = [-1 1 -1 1]; %S
Q(7,[6:8]) = [1 -1 -1]; %SCA
Q(8,[8 22]) = [1 1];  %F

Q(9,[9 10 17 18]) = [-1 1 -1 1];  %QI
Q(10,[11:13]) = [1 -1 1];  %Q
Q(11,[12 13]) = [1 -1];  %QA
Q(12,[9:11]) = [-1 1 1];  %CAQI

Q(13,[14:22]) = [-1 1 1 -1 1 1 -1 1 1];  %N
Q(14,[14:16]) = [1 -1 -1];  %NAI
Q(15,[17:19]) = [1 -1 -1];  %NQI
Q(16,[20:22]) = [1 -1 -1];  %SN

% define simulation parameters
t_max = 1000;
kcat = 0.1; % events per sec per nM
ts = 100; % relative timesale of catalysis and binding/unbinding rates in system
pn_rat = 100; % positive-to-negative activity ratio
ns_rat = 1/100;
% define propensity vector
p_vec = repelem(kcat,size(Q,2))';
p_vec = p_vec.*[ts ts ts ts 1 ts ts 1 ts ts ns_rat ts ts ts ts 1  ts ts 1 ts ts 1]';

% define initial conditions
S0 = 200;
A0 = 1e-8;
AI0 = 0;
C0 = 200;%A0 + AI0;
N0 = A0;
QI0 = AI0/ns_rat/2;

x0_neg = [AI0 0 0 0 C0 S0 0 0 QI0 0 0 0 N0 0 0 0 ];
x0_pos = [AI0 A0 0 0 C0 S0 0 0 QI0 0 0 0 N0 0 0 0 ];
% x0_vec_high = [AI0 pn_rat*A0 0 0 C0 S0 0 0];

% try solver
[t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_neg);
[t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_pos);

f_low_norm = y_low(:,end)/S0*100;
f_high_norm = y_high(:,end)/S0*100;

% Make figure
cmap1 = brewermap(9,'Set2');

ncr_fig = figure;
hold on
plot(t_low,f_low_norm,'Color',cmap1(2,:),'LineWidth',1.5')
plot(t_high,f_high_norm,'Color',cmap1(3,:),'LineWidth',1.5')

plot(t_high,repelem(20,numel(t_high)),'--','Color','black','LineWidth',1)

xlim([100 800])
ylim([0 110])
grid on
xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend('negative sample','positive sample','detection threshold','Location','southeast')
% h = colorbar;
% ylabel(h,'initial target concentration (nM)')
set(gca,'Fontsize',14)
saveas(ncr_fig, [FigPath 'naive_ncr_detection_challenge.png'])
saveas(ncr_fig, [FigPath 'naive_ncr_detection_challenge.pdf'])



% %% test logistic growth approximation
% f_fun = @(params) params(1) / (1 + (params(1) - params(2)) / params(2) * exp(-params(3) * params(4)));
% 
% kon = 1;
% koff = 1;
% kcat = .1;
% 
% r = kon / (kon + koff) * kc;
% t_vec = 0:t_max;
% f_vec_log = NaN(size(t_vec)) ;
% for t = t_vec
%     f_vec_log(t+1) = f_fun([S0 1 r t]);
% end