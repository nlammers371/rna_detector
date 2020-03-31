function dydt = ncr_damper_v1(t,y)

    kon1 = 1;
    koff1 = 1;
    kon2 = 0.1;
    koff2 = 0.1;
    kc = 0.01;
    

    %%%%%%%% Define system of ODEs
    % reactant key:
    %   species 1: QI (caged quencher)
    %   species 2: Q (free quencher)
    %   species 3: QA (quenched activator)
    %   species 4: AI (caged activator)
    %   species 5: A (free activator)
    %   species 6: CA (Cas13:Activator)
    %   species 7: CAQI (Cas13:Activator:Caged quencher)
    %   species 8: CAAI (Cas13:Activator:Caged activator)
    %   species 9: C (free Cas13)
    
    % reaction key:
    %    reaction 1: Q + A -> QA    %kon1
    %    reaction 2: C + A -> CA    %kon1
    %    reaction 3: CA -> C + A    %koff1
    
    %    reaction 4: CA + QI -> CAQI    %kon2
    %    reaction 5: CA + AI -> CAAI    %kon2
    %    reaction 6: CAQI -> CA + QI    %koff2
    %    reaction 7: CAAI -> CA + AI    %koff2
    
    %    reaction 8: CAQI -> CA + Q    %kc
    %    reaction 9: CAAI -> CA + A    %kc
    
    
    % stoichiometry matrix
    Q = zeros(9,9);    
    Q(1,[4 6]) = [-1 1];
    Q(2,[1 8]) = [-1 1];
    Q(3,1) = 1;
    Q(4,[5 7]) = [-1 1];
    Q(5,[1 2 3 9]) = [-1 -1 1 1];
    Q(6,[2 4 5 6 7 8 9]) = [1 -1 -1 1 1 1 1];
    Q(7,[4 6]) = [1 -1];
    Q(8,[5 7]) = [1 -1];
    Q(9,[2 3]) = [-1 1];
    
    % rate vector
    v = zeros(7,1);
    v(1) = kon1*prod(y(Q(:,1)<0));
    v(2) = kon1*prod(y(Q(:,2)<0));
    v(3) = koff1*prod(y(Q(:,3)<0));
    v(4) = kon2*prod(y(Q(:,4)<0));
    v(5) = kon2*prod(y(Q(:,5)<0));
    v(6) = koff2*prod(y(Q(:,6)<0));
    v(7) = koff2*prod(y(Q(:,7)<0));
    v(8) = kc*prod(y(Q(:,8)<0));
    v(9) = kc*prod(y(Q(:,9)<0));
    
    
    % update vector
    dydt = Q*v;
    

end