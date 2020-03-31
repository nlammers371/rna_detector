function [value, isterminal, direction] = myEvent(TIME, Y)

value      = double(any(Y(5) >= Y(4)));
isterminal = 1;   % Stop the integration
direction  = 0;