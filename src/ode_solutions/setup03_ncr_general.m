% Script to generate a general ODE model of NCR reaction

clear
close all
addpath('../utilities')
DataPath = '../../out/ode_studies/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key (atomic entitities)
%   species 1: T1 (target mRNA 1)
%   species 2: T2 (target mRNA 2)
%   species 3: T2 (target mRNA 2 (caged))
%   species 4: G1 (guide mRNA 1)
%   species 5: G2 (guide mRNA 2)
%   species 6: C1 (Cas13 targeted to mRNA 1)
%   species 7: C2 (Cas13 targeted to mRNA 2)
%   species 8: S (dark reporter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modelling approach:
    %(1) Let's allow all atoms to associate and disassociate
    
    %(2) A subset of higher order compounds involving specific Cas13
    %    recognition steps will be allowed to associate and disassociate
    %    with eachother and with atoms
    
    %(3) Catalyst and substrate entities will be specified, giving rise to
    %    a second class of interactions
    

% initialize basic rate symbols
syms kon koff_ns koff_s kcat_ratio kcat_cas sG1 sG2 sT1 sT2 sT2I sSs

% generate list of raction atoms 
atom_string_cell = {'C1' 'C2' 'G1' 'G2' 'T1' 'T2' 'T2I' 'S' };
n_atoms = numel(atom_string_cell);
atom_catalyst_vec = [1 1 0 0 0 0 0 0];
atom_cat_rate_vec = [kcat_cas kcat_cas 0 0 0 0 0 0];
atom_sub_rate_vec = [0 0 sG1 sG2 sT1 sT2 sT2I sSs];

% stable compound reaction components (assuming all non-specific compound
% abundances are negligible)
cas_cat_cmp_string_cell = {'G1:C1','T1:G1:C1','G2:C2','T2:G2:C2'}; % reverse sense from equivalent ns compounds
n_cas_compounds = numel(cas_cat_cmp_string_cell); 
compound_cat_rate_vec = [kcat_cas  kcat_cas*kcat_ratio kcat_cas kcat_cas*kcat_ratio];

% define specific interactions taht do no involve Cas13. Assume that these
% also exist in reasonably high levels within the reaction
stable_cmp_string_cell = {'G1:T1','G2:T2'};
n_stable_compounds = numel(stable_cmp_string_cell);
% catalytic products
sub_indices = find(atom_sub_rate_vec~=0);
cat_product_string_cell = {};
n_products = numel(sub_indices);
for i = sub_indices
    cat_product_string_cell = [cat_product_string_cell {[atom_string_cell{i} 'Z']}];
end

% construct full list of possible compound products

% initialize vectors
n_atoms_total = n_atoms + n_cas_compounds + n_stable_compounds + n_products;
full_atom_cell = [atom_string_cell cas_cat_cmp_string_cell stable_cmp_string_cell cat_product_string_cell];
full_cat_rate_vec = [atom_cat_rate_vec compound_cat_rate_vec zeros(1,n_stable_compounds) zeros(1,n_products)];
full_cat_bin_vec = logical(full_cat_rate_vec~=0);
full_sub_rate_vec = [atom_sub_rate_vec zeros(1,n_cas_compounds) zeros(1,n_stable_compounds) zeros(1,n_products)];
full_sub_bin_vec = logical(full_sub_rate_vec~=0);
full_compound_string_cell = cell(1,n_atoms_total);
cat_rate_vec_sym = [];
cat_product_ind_vec = [];
cat_cat_ind_vec = [];
component_array = NaN(n_atoms_total, 2);
% iterate through all possible combinations
% allow for unique catalytic rate for every catalyst-substrate combination

iter = 1;
for a1 = 1:n_atoms_total
    for a2 = a1:n_atoms_total % allow for self-interactions                
        % generate compound string
        cp_string = [full_atom_cell{a1} '::' full_atom_cell{a2}];        
        spec_flag = any(strcmp(cp_string,stable_cmp_string_cell));
        if ~spec_flag            
            component_array(iter,:) = [a1 a2];
            % add to list
            full_compound_string_cell{iter} = cp_string;

            % check for catalytic potential
            cat_vec = [full_cat_bin_vec(a1) full_cat_bin_vec(a2)];
            cat_indices = find(cat_vec);
            sub_vec = [full_sub_bin_vec(a1) full_sub_bin_vec(a2)];
            sub_indices = find(sub_vec);
            
            if numel(cat_indices) == 1 && numel(sub_indices) == 1                
                if sub_indices == 1
                    cat_rate_vec_sym = [cat_rate_vec_sym full_cat_rate_vec(a2)*full_sub_rate_vec(a1)];
                    pd_ind = find(strcmp(full_atom_cell,[full_atom_cell{a1} 'Z']));
                    cat_product_ind_vec = [cat_product_ind_vec pd_ind];
                    cat_cat_ind_vec = [cat_cat_ind_vec a2];
                else
                    cat_rate_vec_sym = [cat_rate_vec_sym full_cat_rate_vec(a1)*full_sub_rate_vec(a2)];
                    pd_ind = find(strcmp(full_atom_cell,[full_atom_cell{a2} 'Z']));
                    cat_product_ind_vec = [cat_product_ind_vec pd_ind];
                    cat_cat_ind_vec = [cat_cat_ind_vec a1];
                end                
            else
                cat_rate_vec_sym = [cat_rate_vec_sym 0];
                cat_product_ind_vec = [cat_product_ind_vec NaN];
                cat_cat_ind_vec = [cat_cat_ind_vec NaN];
            end
            % increment
            iter = iter + 1;
        end        
    end
end

on_rate_vec_sym = repelem(kon,numel(cat_rate_vec_sym));
off_rate_vec_sym = repelem(koff_ns,numel(cat_rate_vec_sym));

% incorporate the cas catalytic reactions into reference vectors
full_compound_string_cell = [full_compound_string_cell cas_cat_cmp_string_cell];
on_rate_vec_sym = [on_rate_vec_sym repelem(kon,n_cas_compounds)];
off_rate_vec_sym = [off_rate_vec_sym repelem(koff_s,n_cas_compounds)];
cat_rate_vec_sym = [cat_rate_vec_sym repelem(0,n_cas_compounds)];
cat_product_ind_vec = [cat_product_ind_vec NaN(1,n_cas_compounds)];
cat_cat_ind_vec = [cat_cat_ind_vec NaN(1,n_cas_compounds)];

component_array(iter:iter+n_cas_compounds-1,:) = [1 3 ; 5 9 ; 2 4 ;6 11];

% incorporate specific interactions
full_compound_string_cell = [full_compound_string_cell stable_cmp_string_cell];
on_rate_vec_sym = [on_rate_vec_sym repelem(kon,n_stable_compounds)];
off_rate_vec_sym = [off_rate_vec_sym repelem(koff_s,n_stable_compounds)];
cat_rate_vec_sym = [cat_rate_vec_sym repelem(0,n_stable_compounds)];
cat_product_ind_vec = [cat_product_ind_vec NaN(1,n_stable_compounds)];
cat_cat_ind_vec = [cat_cat_ind_vec NaN(1,n_stable_compounds)];

component_array(iter+n_cas_compounds:iter+n_cas_compounds+n_stable_compounds-1,:) = [3 5 ; 4 6];

%% Generate stoichiometry matrix and full-length first order rate vector

% calculate n reactants and n reactions 
n_full_atoms = numel(full_atom_cell);
n_cas_compounds = numel(full_compound_string_cell);
n_reactants = n_cas_compounds + n_full_atoms; % total number of ractions
n_reactions = 2*n_cas_compounds + sum(~isnan(cat_cat_ind_vec));

% initialize stoichiometry matrix
stoichiometry_matrix = zeros(n_reactants,n_reactions);

% strings indicating raction scheme
reaction_string_cell = cell(1,n_reactions);

% initialize vector objects to store reaction rate info
rate_vec_first_order_sym = sym('NA',[1,n_reactions]);
rate_vec_full_sym = sym('NA',[1,n_reactions]);
catalytic_activity_flag_vec = [];
specific_activity_flag_vec = [];
% initialize rate counter
rate_ind = 1;
% iterate through all possible chemical interactions
for i = 1:numel(full_compound_string_cell)
    % check for special cases
    cat_rate = cat_rate_vec_sym(i);  
    off_rate = off_rate_vec_sym(i);
    % extract compound
    cpnd_str = full_compound_string_cell{i};
    cpnd_ind = n_full_atoms + i;
    % get component indices
    component_indices = component_array(i,:);            
    cp1_str = full_atom_cell{component_indices(1)};
    cp2_str = full_atom_cell{component_indices(2)};
    % add association reaction and update Q
    rate_vec_first_order_sym(rate_ind) = on_rate_vec_sym(i);
        
    stoichiometry_matrix(component_indices, rate_ind) = -1;
    stoichiometry_matrix(cpnd_ind,rate_ind) = 1;
    
    % update reaction string array
    reaction_string_cell{rate_ind} = [cp1_str ' + ' cp2_str ' -> ' cpnd_str];
    
    % update special case vec
    catalytic_activity_flag_vec = [catalytic_activity_flag_vec logical(cat_rate~=0)];
    specific_activity_flag_vec = [specific_activity_flag_vec logical(off_rate~=koff_ns)];
    %increment 
    rate_ind = rate_ind + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% add disassociation reaction    
    rate_vec_first_order_sym(rate_ind) = off_rate;

    stoichiometry_matrix(component_indices, rate_ind) = 1;
    stoichiometry_matrix(cpnd_ind,rate_ind) = -1;
    
    % update reaction string array
    reaction_string_cell{rate_ind} = [cpnd_str ' -> ' cp1_str ' + ' cp2_str];
    
    % update special case vec
    catalytic_activity_flag_vec = [catalytic_activity_flag_vec logical(cat_rate~=0)];
    specific_activity_flag_vec = [specific_activity_flag_vec logical(off_rate~=koff_ns)];
    
    % increment 
    rate_ind = rate_ind + 1;   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% check for catalytic activity 
    

    if cat_rate~=0
        % update rate vec and stoichiometry matrix
        rate_vec_first_order_sym(rate_ind) = cat_rate;
        
        % update stoich
        product_ind = cat_product_ind_vec(i);
        product_str = full_atom_cell{product_ind};
        cat_ind = cat_cat_ind_vec(i);
        cat_str = full_atom_cell{cat_ind};        
        
        stoichiometry_matrix([cat_ind product_ind],rate_ind) = 1;
        stoichiometry_matrix(cpnd_ind,rate_ind) = -1;
        
        % update reaction string array
        reaction_string_cell{rate_ind} = [cpnd_str ' -> ' cat_str ' + ' product_str];
        
        % update special case vec
        catalytic_activity_flag_vec = [catalytic_activity_flag_vec logical(cat_rate~=0)];
        specific_activity_flag_vec = [specific_activity_flag_vec logical(off_rate~=koff_ns)];
        
        % increment      
        rate_ind = rate_ind + 1;            
    end   
end

save([DataPath 'ncr_general_setup.mat'])
