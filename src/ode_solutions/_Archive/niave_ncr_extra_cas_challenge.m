% Script to illustrate the sensitivity issue with naive NCR

% To get amplification, we need lots of Cas13, but having lots of Cas13
% gives rise to significant levels of background activity that "looks"
% exactly like activated Cas13 activty wrpt cutting specificity

clear
close all
addpath('../utilities')

FigPath = '../../fig/preliminary_studies/extra_cas/';
mkdir(FigPath)
DataPath = '../../out/preliminary_studies/extra_cas/';
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
%   species 9: CAI
%   species 10: SC
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

%    reaction 9: C + AI -> CAI    %kon 
%    reaction 10: CAI -> C + AI    %koff
%    reaction 11: CAI -> C + A    %kc_basal

%    reaction 12: S + C -> SC    %kon
%    reaction 13: SC -> S + C      %koff 
%    reaction 14: SC -> C + F     %kc_basal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define stoichiometry matrix
Q = zeros(10,14);    
Q(1,[3 4 9 10]) = [-1 1 -1 1];
Q(2,[1 2 5 11]) = [-1 1 1 1];
Q(3,1:8) = [1 -1 -1 1 1 -1 1 1];
Q(4,[3 4 5]) = [1 -1 -1];
Q(5,[1 2 9:14]) = [-1 1 -1 1 1 -1 1 1];
Q(6,[6 7 12:13]) = [-1 1 -1 1];
Q(7,[6:8]) = [1 -1 -1];  
Q(8,[8 14]) = [1 1];  
Q(9,9:11) = [1 -1 -1];  
Q(10,12:14) = [1 -1 -1];  

% define simulation parameters
t_max = 1e5;
kcat = 0.1; % events per sec per nM
ts = 100; % relative timesale of catalysis and binding/unbinding rates in system
pn_rat = 100; % positive-to-negative activity ratio
b_rat = 1/10000; % basal cleavage rate of cas (relative to kcat)
% define propensity vector
p_vec = repelem(kcat,14)'.*[ts ts ts ts 1 ts ts 1 ts ts b_rat ts ts b_rat]';

% define initial conditions
S0 = 200;
A0_low = 0;
A0_high = 1e-2;
AI0 = S0;
exp_list = -9:2:9;
C0_vec = logspace(-8,8,5)*A0_high;

% solve for each Cas13 concentration
results_struct = struct;
for a = 1:numel(exp_list)  
    C0 = A0_high * 10^exp_list(a);
    x0_vec_pos = [AI0 A0_high 0 0 C0 S0 0 0 0 0];
    x0_vec_neg = [AI0 A0_low 0 0 C0 S0 0 0 0 0];
    results_struct(a).x0_vec_high = x0_vec_pos;
    results_struct(a).x0_vec_low = x0_vec_neg;
    
    % try solver
    [t_low,y_low] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_vec_neg);
    [t_high,y_high] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max],x0_vec_pos);

    f_low_norm = y_low(:,8)/S0*100;
    f_high_norm = y_high(:,8)/S0*100;

    % Make figure
    cmap1 = brewermap(9,'Set2');

    ncr_fig = figure;
    hold on
   
    plot(t_low,f_low_norm,'Color',cmap1(3,:),'LineWidth',1.5')
    plot(t_high,f_high_norm,'--','Color',cmap1(2,:),'LineWidth',1.5')

    plot(t_high,repelem(20,numel(t_high)),'--','Color','black','LineWidth',1)

    xlim([0 t_max])
    ylim([0 110])
    grid on
    xlabel('time (seconds)')
    ylabel('fluorescent signal (percent of max)')
    legend('negative sample','positive sample','detection threshold','Location','southeast')
     set(gca,'XScale','log')
    set(gca,'YScale','log')
    % h = colorbar;
    % ylabel(h,'initial target concentration (nM)')
    set(gca,'Fontsize',14)
    saveas(ncr_fig, [FigPath 'naive_ncr_extra_cas_A0_' num2str(A0_high) '_exp' num2str(exp_list(a)) '.png'])
    saveas(ncr_fig, [FigPath 'naive_ncr_extra_cas_A0_' num2str(A0_high) '_exp' num2str(exp_list(a)) '.pdf'])
end
close all
% save([DataPath 'extra_cas_data.mat'])


%% what hapens if we reduce C?
% solve for each Cas13 concentration
results_struct_c0 = struct;
t_max_new = 3600;
for a = 1:numel(A0_high)
    A0_init = A0_high(a);
    C0_init = A0_init / b_rat;
    
    x0_vec_neg = [AI0 A0_low 0 0 C0_init S0 0 0 0 0];
    x0_vec_pos = [AI0  A0_init 0 0 C0_init S0 0 0 0 0];
    results_struct_C0(a).x0_vec_high = x0_vec_pos;
    results_struct_C0(a).x0_vec_low = x0_vec_neg;
    
    % try solver
    [t_low,y_low] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max_new],x0_vec_neg);
    [t_high,y_high] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q),[0 t_max_new],x0_vec_pos);

    f_low_norm = y_low(:,8)/S0*100;
    f_high_norm = y_high(:,8)/S0*100;

    % Make figure
    cmap1 = brewermap(9,'Set2');

    ncr_fig = figure;
    hold on
    plot(t_low,f_low_norm,'Color',cmap1(3,:),'LineWidth',1.5')
    plot(t_high,f_high_norm,'--','Color',cmap1(2,:),'LineWidth',1.5')

    plot(t_high,repelem(20,numel(t_high)),'--','Color','black','LineWidth',1)

    xlim([0 t_max_new])
    ylim([0 110])
    grid on
    xlabel('time (seconds)')
    ylabel('fluorescent signal (percent of max)')
    legend('negative sample','positive sample','detection threshold','Location','southeast')
    % h = colorbar;
    % ylabel(h,'initial target concentration (nM)')
    set(gca,'Fontsize',14)
    saveas(ncr_fig, [FigPath 'naive_ncr_extra_cas_variable_C0_A0' num2str(round(1e8*A0_init)) '.png'])
    saveas(ncr_fig, [FigPath 'naive_ncr_extra_cas_variable_C0_A0' num2str(round(1e8*A0_init)) '.pdf'])
end
close all
save([DataPath 'extra_cas_data.mat'])


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