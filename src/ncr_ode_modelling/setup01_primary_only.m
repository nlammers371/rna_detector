% Script to generate a "simple" ODE model that describes the dynamics of 
% Cas13 detection in the absence of NCR amplification
clear
close all
addpath('../utilities')
DataPath = '../../out/ncr_ode_modeling/';
mkdir(DataPath)

% add specific suffix to indicate version
project_suffix = '_v1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelling approach:
    %(1) Model only "on-pathway" ractions: those between targeted ligans/
    %    substrates and those with catalytic potential. Allow competition from
    %    non-specific interactions to manifest through decreased effective
    %    kon
    
    %(2) Wherever possible, use mass conservation to reduce number of
    %    equations
    
    %(3) Assume uniform association rate (kon) but allow for specific
    %    interactions to be unique to each target-substrate pair
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reactant key (atomic entitities)
%   species 1: A1 (target mRNA 1)
%   species 2: G1 (guide mRNA 1)
%   species 3: NG1 (cleaved guide mRNA 1)
%   species 4: C13 (Cas13)
%   species 5: S (dark reporter)
%   species 6: F (fluorescent reporter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outstanding questions:       
    % (1) Do we need to add a term for Cas13 degredation?
    % (2) What about activator degredation?
    % (3) Is it wrong to assume the same on rate for all association
    %     processes?
    % (4) What is the porper degradation pathway for C13:G1::A1? 
    %     C13 + G1:A1 or C13:G1 + A1?
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define basic reaction components

% Generate list of relevant ractants
atom_list = {'A1','G1','NG1','C13','S','F'};

% Define specific interactions
spec_compound_list = {'C13:G1','G1:A1','C13:G1::A1'}; % double colon denotes where split occurs during dissociation
spec_off_rate_list = sym({'rcg','rga','rga'}); % disassociation rates for these compounds
syms b kc % initialize symbolic variables 
cat_rates =[b*kc, 0, kc]; % b << 1 is backgrounds to activated Cas13 ratio

% Specify substrate-product pairs
substrate_list = {'S','G1'};
product_list = {'F','NG1'};
propensity_list = {'ps','pg1'};

    
% Generate list of catalyst-substrate compounds. These are not specific
% (no base-pairing), but they do have catalytic potential

cat_indices = find(cat_rates~=0);
cat_compound_list = {}; % track compounds
cat_product_list = {}; % track products that result from catalysis in each case
cat_rate_list = [];
for c = cat_indices 
    catalyst = spec_compound_list{c};
    n_colons = contains(catalyst,'::') + 2;
    for s = 1:numel(substrate_list)
        substrate = substrate_list{s};
        cat_compound_list = [cat_compound_list {[catalyst repelem(':',n_colons) substrate]}];
        cat_product_list = [cat_product_list product_list(s)];
        cat_rate_list = [cat_rate_list cat_rates(c)];
    end
end
cat_off_rates = sym(repelem({'rns'},length(cat_compound_list))); % assuming all nonspecific reaction have same off rate


%%% generate stoichiometry matrix and rate vector

% concatenate reactant lists
full_reactant_list = [atom_list spec_compound_list cat_compound_list]; % all reaction elements
full_compound_list = [spec_compound_list cat_compound_list]; % all compounds
full_off_rate_list = [spec_off_rate_list cat_off_rates]; % off rates for each compounds
full_cat_rate_list = [zeros(size(spec_off_rate_list)) cat_rate_list]; % catalytic rates for each compund
full_on_rate_list = sym(repelem({'k'},length(full_off_rate_list))); % association rates

% get dims for stoichiometry matrix (Q)
n_reactants = length(full_reactant_list);
n_reactions = 2*length(full_compound_list) + length(cat_compound_list);

% initialize stoichiometry matrix and rate vector
rate_vec = []; % symbolic list
reaction_list = cell(1,n_reactions); % string cell array
Q = zeros(n_reactants,n_reactions);

% loop through compounds and populate reaction arrays
delim_list = {':','::',':::'};
rate_index = 1;
for f = 1:length(full_compound_list)
    
    %%%%%%%%%%%
    % determine reactants and obtain relevant Q indices    
    
    compound = full_compound_list{f};
    % split reaction
    delim_ind = find([contains(compound,':') contains(compound,'::') contains(compound,':::')],1,'last');
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

% save workspace
save([DataPath 'primary_only' project_suffix '.mat'])
