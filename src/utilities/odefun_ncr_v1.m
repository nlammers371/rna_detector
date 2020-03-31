function dydt = odefun_ncr_v1(t,y)

    % rate parameters 
    koff1 = 9;%1e-3;
    kon1 = 1;%1e-3;
    kc1 = 1;

    koff2 = 9;%1e-3;
    kon2 = 1;%1e-3;
    kc2 = 1;%1e-3;
    v = 1;
    
%     if t > 1
%         error('awrwa')
%     end
    dydt = zeros(6,1);
    % Define update equations
    dydt(1) =  -kon1*y(1)*y(2) + (koff1+kc1)*y(3);% ...
                                    %-kon2*y(1)*y(5) + (koff2+kc2)*y(6) + v*kc2*y(6); % active complex

    dydt(2) =  -kon1*y(1)*y(2) + koff1*y(3); % S

    dydt(3) =  -dydt(1); %AS

%     dydt(4) = kc1*y(3); % cleaved reporter

%     dydt(5) =  -kon2*y(1)*y(5) + koff2*y(6);

%     dydt(6) =  kon2*y(1)*y(5) - (koff2+kc2)*y(6);
end