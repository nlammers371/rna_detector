function F = eq_ode(y)
    kon = 1;
    koff_s = 1;
    koff_ns = 2;
    C013 = 200;
    S0 = 200;
    A013 = 200;
    F = zeros(size(y));
    
    F(4) = kon*y(1)*y(2) - koff_s*y(4);
    F(5) = kon*y(1)*y(3) - koff_ns*y(5);
    F(3) =  y(3) + y(5) - S0;
    F(1) =  F(4) + F(5) + F(1) - C013;
    F(2) = F(4) + F(2) - A013;
end