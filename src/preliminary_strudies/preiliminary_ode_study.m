clear
close all
% Numerical simulations to simple ODEs to convey general pricniples behind
% NCR
addpath('utilities')
FigPath = '../fig/preliminary_studies/';
mkdir(FigPath)

DataPath = '../out/preliminary_studies/';
mkdir(DataPath)


%% NCR
total_time = 3600;

% initial concentration of active complex 
A0_vec = [1e-8 1e-6 1e-4 1e-2 1e0]*10^-9; 

% define valency parameter
valency_vec = [0, 1];

% other initial concentrations (in molar)
S0 = 200e-9; % dark reporter
P0 = 0; % cleaved (fluorescent) reporter
B0 = 200e-9; % caged amplifier
AB0 = 0; % cage-active complex
AS0 = 0;
% rate parameters 
koff1 = 9;%1e-3;
kon1 = 1e7;%1e-3;
kc1 = 1;

koff2 = 9;%1e-3;
kon2 = 1e7;%1e-3;
kc2 = 1;%1e-3;

dt = 0.01;%1/max([kon1 kon2]);
n_steps = round(total_time/dt);
%%
% Define update equations
eqA_ncr = @(c_vec,v) -kon1*c_vec(1)*c_vec(2) + (koff1+kc1)*c_vec(3) ...
                                -kon2*c_vec(1)*c_vec(5) + (koff2+kc2)*c_vec(6) + v*kc2*c_vec(6); % active complex

eqS_ncr = @(c_vec) -kon1*c_vec(1)*c_vec(2) + koff1*c_vec(3); % dark reporter

eqAS_ncr = @(c_vec) kon1*c_vec(1)*c_vec(2) - (koff1+kc1)*c_vec(3);

eqP_ncr = @(c_vec) kc1*c_vec(3); % cleaved reporter

eqB_ncr = @(c_vec) -kon2*c_vec(1)*c_vec(5) + koff2*c_vec(6);

eqAB_ncr = @(c_vec) kon2*c_vec(1)*c_vec(5) - (koff2+kc2)*c_vec(6);


% initialize structure
results_struct = struct;

for v = 1:numel(valency_vec)
    valency = valency_vec(v);
    % initialize vectors to store results
    c_array_ncr = NaN(n_steps,6,numel(A0_vec));
    % loop over starting concentration values
    tic
    for a = 1:numel(A0_vec)
        A0 = A0_vec(a);
        % run simulation
        c_array_ncr(1,:,a) = [A0 S0 AS0 P0 B0 AB0];
        % iterate through time steps
        for n = 2:n_steps
            c_prev = c_array_ncr(n-1,:,a);
            % update
            c_array_ncr(n,1,a) = c_array_ncr(n-1,1,a) + dt*eqA_ncr(c_prev,valency);
            c_array_ncr(n,2,a) = c_array_ncr(n-1,2,a) + dt*eqS_ncr(c_prev);
            c_array_ncr(n,3,a) = c_array_ncr(n-1,3,a) + dt*eqAS_ncr(c_prev);
            c_array_ncr(n,4,a) = c_array_ncr(n-1,4,a) + dt*eqP_ncr(c_prev);
            c_array_ncr(n,5,a) = c_array_ncr(n-1,5,a) + dt*eqB_ncr(c_prev);
            c_array_ncr(n,6,a) = c_array_ncr(n-1,6,a) + dt*eqAB_ncr(c_prev);
        end
    end
    results_struct(v).c_array_ncr = c_array_ncr;
    results_struct(v).valency = valency;    
end
toc

%% Make plots
scale_cell = {'Linear','Log'};
close all
cmap1 = brewermap(numel(A0_vec),'Reds');
time_vec = (1:n_steps)*dt / 60;
ncr_fig = figure;
hold on

for s = 1:2
    for a = 1:numel(A0_vec)
        plot(time_vec,results_struct(2).c_array_ncr(:,4,a)/S0*100,'Color',cmap1(a,:),'LineWidth',1.5')
    end
    xlim([0 5])
    ylim([0 110])
    xlabel('time (minutes)')
    ylabel('fluorescent signal (percent of max)')
    % h = colorbar;
    % ylabel(h,'initial target concentration (nM)')
    set(gca,'Fontsize',14,'YScale',scale_cell{s})
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_v1_' scale_cell{s} '.png'])
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_v1_' scale_cell{s} '.pdf'])
end
%%
cmap2 = brewermap(numel(A0_vec),'Blues');

for s = 1:2
    ncr_fig = figure;
    hold on

    for a = 1:numel(A0_vec)
        plot(time_vec,results_struct(1).c_array_ncr(:,4,a)/S0*100,'Color',cmap2(a,:),'LineWidth',1.5')
    end
    % xlim([0 1e5])
    ylim([0 110])
    xlabel('time (minutes)')
    ylabel('fluorescent signal (percent of max)')
    % h = colorbar;
    % ylabel(h,'initial target concentration (nM)')
    set(gca,'Fontsize',14,'YScale',scale_cell{s})
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_v0_' scale_cell{s} '.png'])
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_v0_' scale_cell{s} '.pdf'])
end

%% Comparison 
% time_vec = 1:n_steps;
a_index = 1;
for s = 1:2
    ncr_fig = figure;
    hold on

    plot(time_vec,results_struct(2).c_array_ncr(:,4,a_index)/S0*100,'Color',cmap1(3,:),'LineWidth',1.5')
    plot(time_vec,results_struct(1).c_array_ncr(:,4,a_index)/S0*100,'Color',cmap2(3,:),'LineWidth',1.5')

    xlim([0 5])
    ylim([0 110])
    xlabel('time (minutes)')
    ylabel('fluorescent signal (percent of max)')
    legend('NCR','Sherlock','Location','northwest')
    % grid on
    % h = colorbar;
    % ylabel(h,'initial target concentration (nM)')
    set(gca,'Fontsize',14,'YScale',scale_cell{s})
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_vs_sher_10aM_' scale_cell{s} '.png'])
    saveas(ncr_fig, [FigPath 'ode_sim_ncr_vs_sher_10aM_' scale_cell{s} '.pdf'])
end

% %% (1) Standard Sherlock
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Define simulation parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % simulation steps
% n_steps = 1e6; 
% 
% % initial concentration (in nanomolar)
% A0 = 1; 
% S0 = 200;
% AS0 = 0;
% P0 = 0;
% 
% % rate parameters 
% koff1 = 1e-3;
% kon1 = 1e-3;
% kc1 = 1e-3;
% 
% % Define update equations
% eqA = @(c_vec) -kon1*c_vec(1)*c_vec(2) + (koff1+kc1)*c_vec(3); % active complex
% 
% eqS = @(c_vec) -kon1*c_vec(1)*c_vec(2) + koff1*c_vec(3); % dark reporter
% 
% eqAS = @(c_vec) kon1*c_vec(1)*c_vec(2) - (koff1+kc1)*c_vec(3);
% 
% eqP = @(c_vec) kc1*c_vec(3); % cleaved reporter
% 
% % initialize vectors to store results
% c_array = NaN(n_steps,4);
% 
% % run simulation
% c_array(1,:) = [A0 S0 AS0 P0];
% 
% for n = 2:n_steps
%     c_prev = c_array(n-1,:);
%     c_array(n,1) = c_array(n-1,1) + eqA(c_prev);
%     c_array(n,2) = c_array(n-1,2) + eqS(c_prev);
%     c_array(n,3) = c_array(n-1,3) + eqAS(c_prev);
%     c_array(n,4) = c_array(n-1,4) + eqP(c_prev);
% end
