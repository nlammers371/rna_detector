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
% reactant key (first order)
%   species 1: A13I (caged activator)
%   species 2: A13 (free activator and/or virus)
%   species 3: C13 (free Cas13)
%   species 4: S (dark reporter)
%   species 5: F (Cleaved (fluorescent) reporter)
%   species 6: N (background nuclease contaminant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate list of raction atoms 
atom_cell = {'C13' 'A13' 'A13I' 'N' 'S' 'F'};
atom_list = sym(atom_cell);
n_atoms = numel(atom_list);

% specify specific compund reactants that will be relatively stable
spec_reactants = {'C13_A13'};

% initialize rates
rate_vec = sym({'kon', 'koff_ns', 'koff_s', 'kcat_high', 'kcat_low'});
syms kon koff_ns koff_s kcat_high kcat_low
% specify reactants with catalytic capacity
% cat_reactants = {'N','C13_A13','C13'};
cat_reactants = {'N','C13'};
% cat_rates = [kcat_high kcat_high kcat_low];  
cat_rates = [kcat_high kcat_low];  
% specify reactants that can be cleaved
cat_targets = {'A13I','S'};
cat_products = {'A13','F'};
% add all possible first order compound products to list
base_reactant_list = atom_list;
order_list = ones(size(atom_list));
base_reactant_cell = atom_cell;
% for i = 1:n_atoms
%     for j = i+1:n_atoms
%         new_cp = [char(atom_list(i)) '_' char(atom_list(j))];    
%         % add to list        
%         base_reactant_list =  [base_reactant_list sym(new_cp)];        
%         base_reactant_cell =  [base_reactant_cell{:} {new_cp}];        
%         order_list = [order_list 2];
%     end
% end

% initialize stoichiometry matrix 
n_basic_reactants = numel(base_reactant_list);
n_reactants = n_basic_reactants + n_basic_reactants * (n_basic_reactants-1) / 2; % total number of ractions
n_reactions = n_basic_reactants * (n_basic_reactants-1) + numel(cat_reactants)*numel(cat_targets); % total number of binding/unbinding reactions
Q_mat = zeros(n_reactants,n_reactions);
rate_vec_sym = sym('NA',[1,n_reactions]);
rate_vec_first_order = sym('NA',[1,n_reactions]);
rate_cell_fun = {};
% now we're ready to iterate through all possible reactions 
full_reactant_list = base_reactant_list;
full_reactant_cell = base_reactant_cell;
rate_ind = 1;
q_ind = 1+n_basic_reactants;
cat_rate_indices = [];
cat_cp_indices = [];

% first do association and disassociation reactions
for i = 1:n_basic_reactants
    for j = i+1:n_basic_reactants
        reactant1 = char(base_reactant_list(i));
        reactant2 = char(base_reactant_list(j));
        new_cp = [reactant1  '_' reactant2];   
        % check to see if product has already been added        
        new_flag = ~(any(strcmp(full_reactant_cell,[reactant1  '_' reactant2])) || ...
                any(strcmp(full_reactant_cell,[reactant2  '_' reactant1])));
        if new_flag
            % add to list        
            full_reactant_list =  [full_reactant_list sym(new_cp)];
            full_reactant_cell =  [full_reactant_cell{:} {new_cp}];
            order_list = [order_list order_list(i)+order_list(j)];
            species_ind = q_ind;
        else
            ordering = find([any(strcmp(full_reactant_cell,[reactant1  '_' reactant2])) ...
                any(strcmp(full_reactant_cell,[reactant2  '_' reactant1]))]);
            if ordering == 1
                species_ind = find(strcmp(full_reactant_cell,[reactant1  '_' reactant2]));
            elseif ordering == 2
                species_ind = find(strcmp(full_reactant_cell,[reactant2  '_' reactant1]));
            else
                error('inconsistent logic')
            end
        end
                
        % add association reaction and update Q
        rate_vec_sym(rate_ind) = base_reactant_list(i)*base_reactant_list(j)*kon;
        rate_vec_fun{rate_ind} = @(y,r) y(i)*y(j)*kon;
        rate_vec_first_order(rate_ind) = kon;
        Q_mat([i j],rate_ind) = -1;
        Q_mat(species_ind,rate_ind) = 1;
        
        % increment 
        rate_ind = rate_ind + 1;
        
        % add disassociation reaction
        if any(ismember(spec_reactants,new_cp))
            rate_vec_sym(rate_ind) = full_reactant_list(species_ind)*koff_s;
            rate_vec_fun{rate_ind} = @(y) y(species_ind) * koff_s;
            rate_vec_first_order(rate_ind) = koff_s;
        else
            rate_vec_sym(rate_ind) = full_reactant_list(species_ind)*koff_ns;
            rate_vec_fun{rate_ind} = @(y) y(species_ind) * koff_ns;
            rate_vec_first_order(rate_ind) = koff_ns;
        end
        Q_mat([i j],rate_ind) = 1;
        Q_mat(species_ind,rate_ind) = -1;
        
        % increment 
        rate_ind = rate_ind + 1;
        
        % check for catalyitic activity
        cat_reactant_flags = strcmp(reactant1,cat_reactants)|strcmp(reactant2,cat_reactants);
        cat_target_flags = strcmp(reactant1,cat_targets)|strcmp(reactant2,cat_targets);
        
        if any(cat_reactant_flags) && any(cat_target_flags)         
            cat_r = cat_reactants(cat_reactant_flags);
            target_r = cat_targets(cat_target_flags);
            product_r = cat_products(cat_target_flags);
            cat_rate = cat_rates(cat_reactant_flags);
            % update rate vec
            rate_vec_sym(rate_ind) = full_reactant_list(species_ind)*cat_rate;
            rate_vec_fun{rate_ind} = @(y) y(species_ind) * cat_rate;
            rate_vec_first_order(rate_ind) = cat_rate;
            % update Q
            product_ind = find(strcmp(full_reactant_cell,product_r)); % catalyzed product
            cat_ind = find(strcmp(full_reactant_cell,cat_r)); % catalyst           
            Q_mat([cat_ind product_ind],rate_ind) = 1;
            Q_mat(species_ind,rate_ind) = -1;
            
            % increment
            cat_rate_indices = [cat_rate_indices rate_ind];
            cat_cp_indices = [cat_cp_indices species_ind];
            rate_ind = rate_ind + 1;            
        end   
        if new_flag
            q_ind = q_ind + 1;
        end
    end
end


%% call ODE solver
t_max = 1e5;
rate_vec_val = [1 1e4 1 1 1];
p_vec = NaN(size(rate_vec_first_order));
% generate valued long rate vec
for r = 1:numel(rate_vec_val)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);
%     rep_indices = rate_vec_first_order==rep_rate;
    p_vec(ri) = double(subs(rr,rr,rate_vec_val(r)));       
end

atom_cell = {'C13' 'A13' 'A13I' 'N' 'S' 'F'};
y0_vec = zeros(1,size(Q_mat,1));
y0_vec(1:5) = [100 1e-8 0 0 200];


[t ,y] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q_mat),[0 t_max],y0_vec);