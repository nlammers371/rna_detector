function dydt = ncr_naive(t,y)

    koff_rs = 1;%param_struct.koff1;%1e-3;
    kon_rs = 1;%param_struct.kon1;%1e-3;
    kc_rs = 1;%param_struct.kc1;
    

    %%%%%%%% Define system of ODEs
    % key:
    %    (reaction 1) S + R -> SR
    %    (reaction 2) SR -> R+S
    %    (reaction 3) SR -> R+F
    
    % stoichiometry matrix
    Q = zeros(4,3);    
    Q(1:2,1) = -1;
    Q(3,1) = 1;
    Q(1:2,2) = 1;
    Q(3,2) = -1;
    Q([1 4],3) = 1;
    Q(3,3) = -1;
    
    % rate vector
    v = zeros(3,1);
    v(1) = kon_rs*y(1)*y(2);
    v(2) = koff_rs*y(3);
    v(3) = kc_rs*y(3);
    
    % update vector
    dydt = Q*v;
    
%     % Class 1: reporter substrate binding/unbinding
%     dydt(1) = -S*(kon_rs*R + kon_ns*N + kon_as*A + kon_cs*C) + ...
%                     koff_rs*RS + koff_ns*NS + koff_as*AS + koff_cs*CS; % free substrate
%     dydt(2) = kon_rs*R*S - koff_rs*dydt(2)
%     % (2) cage cleavage
%     
%     
%     
%     dydt = zeros(6,1);
%     % Define update equations
%     dydt(1) =  -kon1*y(1)*y(2) + (koff1+kc1)*y(3);% ...
%                                     %-kon2*y(1)*y(5) + (koff2+kc2)*y(6) + v*kc2*y(6); % active complex
% 
%     dydt(2) =  -kon1*y(1)*y(2) + koff1*y(3); % S
% 
%     dydt(3) =  -dydt(1); %AS

%     dydt(4) = kc1*y(3); % cleaved reporter

%     dydt(5) =  -kon2*y(1)*y(5) + koff2*y(6);

%     dydt(6) =  kon2*y(1)*y(5) - (koff2+kc2)*y(6);
end