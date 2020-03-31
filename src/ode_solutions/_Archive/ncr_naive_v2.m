function dydt = ncr_naive_v3(t,y)

    kon = 1;
    koff = 1;
    kc = 1;
    

    %%%%%%%% Define system of ODEs
    % reactant key:   
    %   species 1: AI (caged activator)
    %   species 2: A (free activator)
    %   species 3: CA (Cas13:Activator)    
    %   species 4: CAAI (Cas13:Activator:Caged activator)
    %   species 5: C (free Cas13)
    %   species 6: S (dark reporter)
    %   species 7: SCA (dark reporter  + Cas13:activator)
    %   species 8: F (fluorescent reporter)
    
    % reaction key:    
    %    reaction 1: C + A -> CA    %kon1
    %    reaction 2: CA -> C + A    %koff1
        
    %    reaction 3: CA + AI -> CAAI    %kon2    
    
    %    reaction 4: CAAI -> CA + AI    %koff2        
    %    reaction 5: CAAI -> CA + A    %kc2
    
    %    reaction 6: S + CA -> SCA    %kon3
    %    reaction 7: SCA -> S + CA      %koff3  
    %    reaction 8: SCA -> CA + F     %kc3      
    
    % stoichiometry matrix
    Q = zeros(8,8);    
    Q(1,[3 4]) = [-1 1 ];
    Q(2,[1 2 5]) = [-1 1 1];
    Q(3,1:8) = [1 -1 -1 1 1 -1 1 1];
    Q(4,[3 4 5]) = [1 -1 -1];
    Q(5,[1 2]) = [-1 1];
    Q(6,[6 7]) = [-1 1];
    Q(7,[6:8]) = [1 -1 -1];  
    Q(8,8) = [1];  
    
    % porpensity vector
    p_vec = [kon koff kon koff kc kon koff kc]';
    % rate vector
    v_vec = zeros(size(Q,2),1);
    % apply mass action
    for i = 1:numel(v_vec)
        v_vec(i) = p_vec(i)*prod(y(Q(:,i)<0));
    end
%     if t > 96
%         error('sigh')
%     end
    % update vector
    dydt = Q*v_vec;
    

end