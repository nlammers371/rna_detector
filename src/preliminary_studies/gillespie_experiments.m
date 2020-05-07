clear
close all

% set experimental time resolution and total time 
T = 360; % seconds
% reaction rates (events per times tep)
kon = 1;
koff = 1;
kp = 0.1;
% definte total concentrations of substrate and enzyme
ET = 1e4;
ST = 1e4;
% define initial vectors for substrate and enzyme
e0 = ET;
s0 = ST;
P0 = 0;
ES = 0;

s_n_vec = [ST];
e_n_vec = [ET];
p_n_vec = [0];
jump_time_vec = [0];

% simulate
t_curr = 0;
tic
while t_curr < T
    % binding reactions first
    kon_curr = s_n_vec(end)*e_n_vec(end)*kon;
    koff_curr = (ET-e_n_vec(end))*koff;
    kp_curr = (ET-e_n_vec(end))*kp;
    tau = 1/(kon_curr+koff_curr+kp_curr); 
    if isinf(tau)
        break
    end
    % sample next time step and state
    next_time = exprnd(tau);
    next_reaction = randsample(1:3,1,true,[kon_curr koff_curr kp_curr]);
    t_curr = t_curr + next_time;
    
    if t_curr < T
        jump_time_vec = [jump_time_vec t_curr];
        
        if next_reaction == 1
            e_n_vec = [e_n_vec e_n_vec(end)-1];
            s_n_vec = [s_n_vec s_n_vec(end)-1];
            p_n_vec = [p_n_vec p_n_vec(end)];
        elseif next_reaction == 2
            e_n_vec = [e_n_vec e_n_vec(end)+1];
            s_n_vec = [s_n_vec s_n_vec(end)+1];
            p_n_vec = [p_n_vec p_n_vec(end)];
        elseif next_reaction == 3
            e_n_vec = [e_n_vec e_n_vec(end)+1];
            p_n_vec = [p_n_vec p_n_vec(end)+1];
            s_n_vec = [s_n_vec s_n_vec(end)];
        end
    end
end
toc

%% Tau leaping method
epsilon = 0.03;
dt_init = .01/kon/ET;
n_steps = round(1 / dt);

% initialize array
reaction_array = NaN(3,n_steps);
reaction_array(:,1) = [ST ; ET ; 0];
% encode stoichiometry of each reaction
stoich_array = [-1, 1, 0 ; 
                -1, 1, 1 ;
                 0, 0, 1];


for n = 2:500
    % enzyme binds to substrate 
    kon_curr = reaction_array(1,n-1)*reaction_array(2,n-1)*kon;    
    % enzyme unbinds from substrate
    koff_curr = (ET-reaction_array(2,n-1))*koff;    
    % cleavage
    kp_curr = (ET-reaction_array(2,n-1))*kp;
    
    % get expected N events
    n_on = poissrnd(kon_curr*dt);
    n_kp = poissrnd(kp_curr*dt);
    n_off = poissrnd(koff_curr*dt);
    % update
    update = [n_on, n_off, n_kp] .* stoich_array;
    reaction_array(:,n) = reaction_array(:,n-1) + sum(update,2);
    
end