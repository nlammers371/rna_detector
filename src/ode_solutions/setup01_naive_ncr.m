% Script to illustrate the sensitivity issue with naive NCR

% Even "negative" samples with low level background activity will trigger
% the cascade, and the timing with which the signals for positve and
% negative samples is quite close, even given relatively large initial
% differences

clear
close all
addpath('../utilities')
DataPath = '../../out/ode_studies/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key (first order)
%   species 1: A13I (caged activator)
%   species 2: A13 (free activator and/or virus)
%   species 3: C13 (free Cas13)
%   species 4: S (dark reporter)
%   species 5: F (Cleaved (fluorescent) reporter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize rates
rate_vec = sym({'kon', 'koff_ns', 'koff_s', 'kcat_high', 'kcat_low'});
syms kon koff_ns koff_s kcat_high kcat_low

% generate list of raction atoms 
atom_cell = {'C13' 'A13' 'A13I' 'S' 'F'};
atom_list = sym(atom_cell);
n_atoms = numel(atom_list);

% compound reaction components
compound_cell = {'C13_A13','C13_A13I','C13_S','C13_A13__A13I','C13_A13__S'};
off_rate_vec = [koff_s koff_ns koff_ns koff_ns koff_ns];
compound_list = sym(compound_cell);
n_compounds = numel(compound_cell); 

% define compounds with catalytic activity
cat_compounds =  {'C13_A13I','C13_S','C13_A13__A13I','C13_A13__S'};
cat_rates = [kcat_low kcat_low kcat_high kcat_high];
cat_products = {'A13','F','A13','F'};

% initialize stoichiometry matrix 
full_reactant_cell = [atom_cell compound_cell];
full_reactant_list = [atom_list compound_list];
n_reactants = n_atoms + n_compounds; % total number of ractions
n_reactions = 2*n_compounds + numel(cat_compounds);
Q_mat = zeros(n_reactants,n_reactions);
         
rate_vec_first_order = sym('NA',[1,n_reactions]);
rate_cell_fun = {};
cat_rate_indices = [];
cat_cp_indices = [];
rate_ind = 1;
% first do association and disassociation reactions
for i = 1:numel(compound_cell)
    % extract compound
    cp_str = compound_cell{i};
    cp_ind = i + n_atoms;
    % extract reactant components
    cp3_flag = contains(cp_str,'__');
    if cp3_flag
        reactants = strsplit(cp_str,'__');
    else
        reactants = strsplit(cp_str,'_');
    end
    r_ind1 = find(strcmp(full_reactant_cell,reactants(1)));
    r_ind2 = find(strcmp(full_reactant_cell,reactants(2)));
    

    % add association reaction and update Q
    rate_vec_sym(rate_ind) = full_reactant_list(r_ind1)*full_reactant_list(r_ind2)*kon;
    rate_vec_fun{rate_ind} = @(y,r) y(r_ind1)*y(r_ind2)*kon;
    rate_vec_first_order(rate_ind) = kon;
    Q_mat([r_ind1 r_ind2],rate_ind) = -1;
    Q_mat(cp_ind,rate_ind) = 1;

    
    %%%%%%%%%%% increment 
    rate_ind = rate_ind + 1;

    % add disassociation reaction
    off_rate = off_rate_vec(i);
    
    rate_vec_sym(rate_ind) = full_reactant_list(cp_ind)*off_rate;
    rate_vec_fun{rate_ind} = @(y) y(cp_ind) * off_rate;
    rate_vec_first_order(rate_ind) = off_rate;

    Q_mat([r_ind1 r_ind2],rate_ind) = 1;
    Q_mat(cp_ind,rate_ind) = -1;

    % increment 
    rate_ind = rate_ind + 1;   
    
    % check for catalytic activity 
    cat_ind = find(strcmp(cat_compounds,cp_str));
    if ~isempty(cat_ind)

        % update rate vec
        rate_vec_sym(rate_ind) = full_reactant_list(cp_ind)*cat_rates(cat_ind);
        rate_vec_fun{rate_ind} = @(y) y(cp_ind) * cat_rates(cat_ind);
        rate_vec_first_order(rate_ind) = cat_rates(cat_ind);
        % update Q
        product_str = cat_products{cat_ind};
        product_ind = find(strcmp(full_reactant_cell,product_str)); % catalyzed product
        
        Q_mat([r_ind1 product_ind],rate_ind) = 1;
        Q_mat(cp_ind,rate_ind) = -1;

        % increment
        cat_rate_indices = [cat_rate_indices rate_ind];
        cat_cp_indices = [cat_cp_indices cp_ind];
        rate_ind = rate_ind + 1;            
    end   
end

save([DataPath 'naive_ncr_setup.mat'])

%% (1) make initial plots to illustrate the promise and challenge




