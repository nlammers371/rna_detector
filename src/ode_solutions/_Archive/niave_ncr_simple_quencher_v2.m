% Script to explore viability of simple quencher model
% Assuming Michaelis emnton to simplify things
clear
close all
addpath('../utilities')

FigPath = '../fig/ode_studies/quencher/';
mkdir(FigPath)
DataPath = '../out/ode_studies/quencher/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key:   
%   species 1: A (activator)
%   species 2: AI (caged activator)
%   species 3: C (Cas13)    
%   species 4: N (Background shit)
%   species 5: Q (Quencher)
%   species 6: QI (caged quencher)
%   species 7: S (dark reporter)
%   species 8: F (Cleaved (fluorescent) reporter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reaction key:    
%    reaction 1: AI -> A    %kon   
%    reaction 2: QI -> Q    %kon   
%    reaction 3:  S -> F    %kon   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Kd = 0.5;

% define stoichiometry matrix
Q = zeros(8,3);    

Q(1,1) = 1; % A
Q(2,1) = -1; %AI
Q(5,2) = 1; % Q
Q(6,2) = -1; % QI
Q(7,3) = -1; %C
Q(8,3) = 1; %S

% define simulation parameters
t_max = 1e5;
kcat = 0.1; % events per sec per nM
ts = 100; % relative timesale of catalysis and binding/unbinding rates in system
pn_rat = 100; % positive-to-negative activity ratio
ns_rat = 1/10000;


% define initial conditions
S0 = 200;
A0 = 1e-4;
AI0 = 200;
C0 = 200;%A0 + AI0;
N0 = A0;
Q0 = S0*1e2;
QI0 = 0;%AI0/ns_rat;

x0_neg = [0 AI0 C0 A0+N0 Q0 QI0 S0 0];
x0_pos = [A0 AI0 C0 N0 0 QI0 S0 0];
% x0_vec_high = [AI0 pn_rat*A0 0 0 C0 S0 0 0];

% try solver
% [t_neg,y_neg] = ode15s(@(t,y) ncr_solver_simp(t,y,Kd,kcat,ns_rat,Q),[0 t_max],x0_neg);
[t_pos,y_pos] = ode15s(@(t,y) ncr_solver_simp(t,y,Kd,kcat,ns_rat,Q),[0 t_max],x0_pos);

f_neg_norm = y_neg(:,end)/S0*100;
f_pos_norm = y_pos(:,end)/S0*100;

% Make figure
cmap1 = brewermap(9,'Set2');

ncr_fig = figure;
hold on
plot(t_neg/t_max,f_neg_norm,'Color',cmap1(3,:),'LineWidth',1.5')
plot(t_pos/t_max,f_pos_norm,'--','Color',cmap1(2,:),'LineWidth',1.5')

plot(t_pos/t_max,repelem(20,numel(t_pos)),'--','Color','black','LineWidth',1)

xlim([0 1])
ylim([0 110])
grid on
xlabel('time (au)')
set(gca,'XScale','log')
ylabel('fluorescent signal (percent of max)')
legend('negative sample','positive sample','detection threshold','Location','southwest')
% h = colorbar;
% ylabel(h,'initial target concentration (nM)')
set(gca,'Fontsize',14)
saveas(ncr_fig, [FigPath 'naive_ncr_quench_Q0' num2str(Q0) '_QI0' num2str(QI0) '.png'])
saveas(ncr_fig, [FigPath 'naive_ncr_quench_Q0' num2str(Q0) '_QI0' num2str(QI0) '.pdf'])



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