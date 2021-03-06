% Script to generate a general ODE model of NCR reaction

clear
close all
addpath('../utilities')
DataPath = '../../out/ncr_ode_modeling/';
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
        
%%%%%%%%%%%%%%%%%   
%%%% (1) Define basic reaction components
%%%%%%%%%%%%%%%%%

% Generate list of relevant ractants
atom_list = {'C13','G1','G2','A1','A2','IA2','I','S','F','N'};

%%%%%%%%%%%%%%%%%
% (1) Define specific interactions
%%%%%%%%%%%%%%%%%

spec_compound_list = {'C13:G1','C13:G2','G1:A1','G2:A2','C13:G1::A1','C13:G2::A2','[C13:G1::A1]','[C13:G2::A2]','I:A2'};
syms k rcga
spec_on_rate_list = [k, k, k, k, k, k, 0, 0, k]; % association rates for these compounds
syms rcg rga kd_cga ria
spec_off_rate_list = [rcg,rcg,rga,rga,k*kd_cga,k*kd_cga,k*kd_cga,k*kd_cga,ria]; % disassociation rates for these compounds
syms b kc
cat_rates =[b*kc, b*kc, 0 ,0, 0, 0, kc, kc, 0];
syms kact
mod_rate_list =[0, 0, 0 ,0, 0, 0, kact, kact, 0];

%%%%%%%%%%%%%%%%%  
%%%% (2) Generate catalyst-substrate compounds
%%%%%%%%%%%%%%%%%

% Specify substrate-product pairs
substrate_list = {'S','G1','G2','IA2'};
product_list = {'F','N','N','I:A2'};
% propensity_list = {'ps','pg1','pg2','pi'};

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



%%%%%%%%%%%%%%%%%
%%% (3) generate stoichiometry matrix and rate vector
%%%%%%%%%%%%%%%%%

% concatenate reactant lists
full_reactant_list = [atom_list spec_compound_list cat_compound_list];
full_compound_list = [spec_compound_list cat_compound_list];
full_off_rate_list = [spec_off_rate_list cat_off_rates];
full_mod_rate_list = [mod_rate_list zeros(size(cat_off_rates))];
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
    
    % check for brackets (is it a modified compound?)
    mod_flag = strcmp(compound(1),'[')&&strcmp(compound(end),']');
    cp_mod = regexprep(compound,'\[|\]','');
    
    % split reaction
    fnd1 = strfind(compound,':');
    fnd2 = strfind(compound,'::');
    fnd3 = strfind(compound,':::');
    delim_ind = find([any(horzcat(fnd1(:))) any(horzcat(fnd2(:))) any(horzcat(fnd3(:)))],1,'last');
    if mod_flag
      reactants = strsplit(cp_mod,delim_list{delim_ind});
    else
      reactants = strsplit(compound,delim_list{delim_ind});
    end
    
    % get indices 
    compound_index = find(strcmp(full_reactant_list,{compound}));
    reactant_ind1 = find(strcmp(full_reactant_list,reactants{1}));
    reactant_ind2 = find(strcmp(full_reactant_list,reactants{2}));
    reactant_indices = [reactant_ind1 reactant_ind2];
    if any([isempty(compound_index) isempty(reactant_ind1) isempty(reactant_ind2)])
        error('problem with reactant indexing')
    end
    
    % check to see if compound is formed through association or
    % modification
    mod_rate = full_mod_rate_list(f);
    
    if mod_rate ~= 0 && ~mod_flag
      error('indonsistent mod labeling')
    elseif mod_rate ~= 0 && mod_flag
      parent_index = find(strcmp(full_reactant_list,{cp_mod}));
      %%%%%%%%%%%
      % modification step
      Q(compound_index,rate_index) = 1;
      Q(parent_index,rate_index) = -1;
      % add reaction string
      reaction_list{rate_index} = [cp_mod ' -> ' compound];
      % add reaction rate
      rate_vec = [rate_vec full_mod_rate_list(f)];
      % increment
      rate_index = rate_index + 1;
    elseif full_on_rate_list(f) ~= 0
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
    else
      error('issue with mod rate assignment')
    end
    
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
save([DataPath 'ncr_v2.mat'])