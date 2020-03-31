% Script to illustrate the sensitivity issue with naive NCR

% Even "negative" samples with low level background activity will trigger
% the cascade, and the timing with which the signals for positve and
% negative samples is quite close, even given relatively large initial
% differences

clear
close all
addpath('../utilities')

FigPath = '../fig/preliminary_studies/';
mkdir(FigPath)
DataPath = '../out/preliminary_studies/';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define stoichiometry matrix
Q = zeros(8,8);    
Q(1,[3 4]) = [-1 1 ];
Q(2,[1 2 5]) = [-1 1 1];
Q(3,1:8) = [1 -1 -1 1 1 -1 1 1];
Q(4,[3 4 5]) = [1 -1 -1];
Q(5,[1 2]) = [-1 1];
Q(6,[6 7]) = [-1 1];
Q(7,[6:8]) = [1 -1 -1];  
Q(8,8) = [1];  

% define simulation parameters
t_max = 1000;
kcat = 0.1; % events per sec per nM
ts = 100; % relative timesale of catalysis and binding/unbinding rates in system
pn_rat = 100; % positive-to-negative activity ratio

% define propensity vector
p_vec = [ts*kcat ts*kcat ts*kcat ts*kcat kcat ts*kcat ts*kcat kcat]';

% define initial conditions
S0 = 200;
A0 = 1e-8;
AI0 = S0;
C0 = A0 + AI0;
x0_vec = [AI0 A0 0 0 C0 S0 0 0];
x0_vec_high = [AI0 pn_rat*A0 0 0 C0 S0 0 0];

% try solver
[t_low,y_low] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_vec);
[t_high,y_high] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_vec_high);

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