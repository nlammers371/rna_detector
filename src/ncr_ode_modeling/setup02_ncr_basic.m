% Script to generate a general ODE model of NCR reaction

clear
close all
addpath('../utilities')
DataPath = '../../out/ode_studies_v2/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key (atomic entitities)
%   species 1: T1 (target mRNA 1)
%   species 2: T2 (target mRNA 2)
%   species 3: IT2 (target mRNA 2 (caged))
%   species 4: G1 (guide mRNA 1)
%   species 5: G2 (guide mRNA 2)
%   species 6: C1 (Cas13 targeted to mRNA 1)
%   species 8: S (dark reporter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modelling approach:
    %(1) Model only "on-pathway" ractions: those between targeted ligans/
    %    substrates and those with catalytic potential. Allow competition from
    %    non-specific interactions to manifest through decreased effective
    %    kon
    
    %(2) Wherever possible, use mass conservation to reduce number of
    %    equations
    
    %(3) Assume uniform association rate (kon) but allow for specific
    %    interactions to be unique to each target-substrate pair
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% (1) Define basic reaction components

% Generate list of relevant ractants
atom_list = {'C13','G1','G2','A1','A2','IA2','I','S','F','N'};

% Define specific interactions
spec_compound_list = {'C13:G1','C13:G2','G1:A1','G2:A2','C13:G1::A1','C13:G2::A2','I:A2'};
syms k rcga
spec_on_rate_list = [k, k, k, k, rcga, rcga, k]; % association rates for these compounds
syms rcg rga kd_cga ria
spec_off_rate_list = [rcg,rcg,rga,rga,rcga*kd_cga, rcga*kd_cga, ria]; % disassociation rates for these compounds
syms b kc
cat_rates =[b*kc, b*kc, 0 ,0, kc, kc, 0];

syms kdeg
unstable_atom_list = {'C13', 'C13:G1', 'C13:G1::A1', 'C13:G1::A1', 'C13:G2', 'C13:G2::A2', 'C13:G2::A2'};
degredation_product = {{'nC13'}, {'nC13', 'G1'}, {'nC13', 'G1:A1'}, {'C13', 'G1:A1'}, {'nC13', 'G2'}, {'nC13', 'G2:A2'}, {'C13', 'G2:A2'}};
deg_rates = {kdeg, kdeg, kdeg, rcga*kd_cga, kdeg, kdeg, rcga*kd_cga};

% Specify substrate-product pairs
substrate_list = {'S','G1','G2','IA2'};
product_list = {'F','N','N','I:A2'};
propensity_list = {'ps','pg1','pg2','pi'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% (2) Generate catalyst-substrate compounds

cat_indices = find(cat_rates~=0);
cat_compound_list = {};
cat_product_list = {};
cat_rate_list = [];
for c = cat_indices 
    catalyst = spec_compound_list{c};
    fnd = strfind(catalyst,'::');
    n_colons = any(horzcat(fnd(:))) + 2;
    for s = 1:numel(substrate_list)
        substrate = substrate_list{s};
        cat_compound_list = [cat_compound_list {[catalyst repelem(':',n_colons) substrate]}];
        cat_product_list = [cat_product_list product_list(s)];
        cat_rate_list = [cat_rate_list cat_rates(c)];
    end
end
cat_off_rates = sym(repelem({'rns'},length(cat_compound_list)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (3) generate stoichiometry matrix and rate vector

% concatenate reactant lists
full_reactant_list = [atom_list spec_compound_list cat_compound_list];
full_compound_list = [spec_compound_list cat_compound_list];
full_off_rate_list = [spec_off_rate_list cat_off_rates];
full_cat_rate_list = [zeros(size(spec_off_rate_list)) cat_rate_list];
full_on_rate_list = [spec_on_rate_list repelem({k},length(full_off_rate_list))];

% get dims for stoichiometry matrix (Q)
n_reactants = length(full_reactant_list);
n_reactions = 2*length(full_compound_list) + length(cat_compound_list);

% initialize stoichiometry matrix and rate vector
rate_vec = [];
reaction_list = {};
Q = zeros(n_reactants,n_reactions);

% loop through compounds and populate reaction arrays
delim_list = {':','::',':::'};
rate_index = 1;
for f = 1:length(full_compound_list)
    %%%%%%%%%%%
    % determine reactants and obtain relevant Q indices
    
    compound = full_compound_list{f};
    % split reaction
    fnd1 = strfind(compound,':');
    fnd2 = strfind(compound,'::');
    fnd3 = strfind(compound,':::');
    delim_ind = find([any(horzcat(fnd1(:))) any(horzcat(fnd2(:))) any(horzcat(fnd3(:)))],1,'last');
    reactants = strsplit(compound,delim_list{delim_ind});
    % get indices 
    compound_index = find(strcmp(full_reactant_list,{compound}));
    reactant_ind1 = find(strcmp(full_reactant_list,reactants{1}));
    reactant_ind2 = find(strcmp(full_reactant_list,reactants{2}));
    reactant_indices = [reactant_ind1 reactant_ind2];
    if any([isempty(compound_index) isempty(reactant_ind1) isempty(reactant_ind2)])
        error('problem with reactant indexing')
    end
    
    %%%%%%%%%%%
    % Binding reaction
    Q(compound_index,rate_index) = 1;
    Q(reactant_indices,rate_index) = -1;
    % add reaction string
    reaction_list{rate_index} = [reactants{1} ' + ' reactants{2} ' -> ' compound];
    % add reaction rate
    rate_vec = [rate_vec full_on_rate_list(f)];
    % increment
    rate_index = rate_index + 1;
    
    %%%%%%%%%%%
    % Unbinding reaction
    Q(compound_index,rate_index) = -1;
    Q(reactant_indices,rate_index) = 1;
    % add reaction string
    reaction_list{rate_index} = [compound ' -> ' reactants{1} ' + ' reactants{2}];
    % add reaction rate
    rate_vec = [rate_vec full_off_rate_list(f)];
    % increment
    rate_index = rate_index + 1;
    
    %%%%%%%%%%%
    % Catalytic reaction (if applicable)
    if full_cat_rate_list(f)~=0
        sub_index = find(strcmp(substrate_list,reactants{2}));
        product_index = find(strcmp(full_reactant_list,product_list(sub_index)));
        % update Q
        Q(compound_index,rate_index) = -1;
        Q([reactant_indices(1) product_index],rate_index) = 1;
        % add reaction string
        reaction_list{rate_index} = [compound ' -> ' reactants{1} ' + ' product_list{sub_index}];
        % add reaction rate
        rate_vec = [rate_vec full_cat_rate_list(f)];
        % increment
        rate_index = rate_index + 1;
    end
end

for f = 1:length(unstable_atom_list)
    compound = unstable_atom_list{f};
    product = degredation_product{f};
    if length(product) == 1
         % get indices 
        compound_index = find(strcmp(full_reactant_list,{compound}));
        product_ind = find(strcmp(full_reactant_list,product{1}));
        %%%%%%%%%%%
        % degredation reaction
        Q(compound_index,rate_index) = -1;
        Q(product_ind,rate_index) = 1;
        % add reaction string
        reaction_list{rate_index} = [compound ' -> ' product{1}];
        % add reaction rate
        rate_vec = [rate_vec deg_rates(f)];
        % increment
        rate_index = rate_index + 1;
    elseif length(product) == 2
             % get indices 
        compound_index = find(strcmp(full_reactant_list,{compound}));
        product_ind1 = find(strcmp(full_reactant_list,product{1}));
        product_ind2 = find(strcmp(full_reactant_list,product{2}));
        product_indices = [product_ind1 product_ind2];
        %%%%%%%%%%%
        % degredation reaction
        Q(compound_index,rate_index) = -1;
        Q(product_indices,rate_index) = 1;
        % add reaction string
        reaction_list{rate_index} = [compound ' -> ' product{1} ' + ' product{2}];
        % add reaction rate
        rate_vec = [rate_vec deg_rates(f)];
        % increment
        rate_index = rate_index + 1;
    end
    
end

reaction_list'
% save workspace
save([DataPath 'ncr_basic.mat'])
