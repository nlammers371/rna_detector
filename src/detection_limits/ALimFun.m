
function A1Lim = ALimFun(S0, B0, RNP1, RNP2, bc, k, wc,rcg)


A1Lim = 8.*(rcg+k.*RNP1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.* ...
  RNP1).^(1/2)).^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.* ...
  RNP2).^(1/2)).*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^( ...
  1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+( ...
  -1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))).*((1/2).*k.^(-1).*(rcg+ ...
  2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).*(1000.*rcg+ ...
  k.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.* ...
  k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2) ...
  .*(rcg+4.*k.*RNP2).^(1/2))+S0))+(-1/2).*(k.^(-2).*(rcg+2.*k.*RNP2+ ...
  (-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).^2.*((-4).*k.^2.*(B0+( ...
  1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1) ...
  .^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+ ...
  4.*k.*RNP2).^(1/2))).*S0+(1000.*rcg+k.*(B0+(1/2).*k.^(-1).*(rcg+ ...
  2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^( ...
  -1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))+ ...
  S0)).^2)).^(1/2)).^(-1).*(1+((1/2).*bc.*k.^(-1).*(B0+(1/2).*k.^( ...
  -1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2))+( ...
  1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2) ...
  .^(1/2))).^(-1).*((1/4).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2) ...
  .*(rcg+4.*k.*RNP2).^(1/2)).^2+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1) ...
  .*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).*(1000.*rcg+k.*(B0+(1/2).* ...
  k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2) ...
  )+S0))+(-1/2).*(k.^(-2).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+ ...
  4.*k.*RNP2).^(1/2)).^2.*((-4).*k.^2.*(B0+(1/2).*k.^(-1).*(rcg+2.* ...
  k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1) ...
  .*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))).*S0+ ...
  (1000.*rcg+k.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2) ...
  .*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).* ...
  rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))+S0)).^2)).^(1/2))+(1/2).*bc.* ...
  k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2) ...
  ).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).^( ...
  -1).*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+ ...
  4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^( ...
  1/2).*(rcg+4.*k.*RNP2).^(1/2))).^(-1).*((1/4).*k.^(-1).*(rcg+2.* ...
  k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2)).*(rcg+2.*k.* ...
  RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))+(1/2).*k.^(-1).*( ...
  rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).*(1000.* ...
  rcg+k.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+ ...
  4.*k.*RNP2).^(1/2))+S0))+(-1/2).*(k.^(-2).*(rcg+2.*k.*RNP2+(-1).* ...
  rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).^2.*((-4).*k.^2.*(B0+(1/2).* ...
  k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2) ...
  )+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.* ...
  RNP2).^(1/2))).*S0+(1000.*rcg+k.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.* ...
  RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*( ...
  rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))+S0)).^2) ...
  ).^(1/2))+2.*(rcg+k.*RNP2).*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+( ...
  -1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.* ...
  k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))).*S0.*((1/2).* ...
  k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2) ...
  ).*(1000.*rcg+k.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^( ...
  1/2).*(rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+( ...
  -1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))+S0))+(-1/2).*(k.^(-2).*( ...
  rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2)).^2.*(( ...
  -4).*k.^2.*(B0+(1/2).*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*( ...
  rcg+4.*k.*RNP1).^(1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).* ...
  rcg.^(1/2).*(rcg+4.*k.*RNP2).^(1/2))).*S0+(1000.*rcg+k.*(B0+(1/2) ...
  .*k.^(-1).*(rcg+2.*k.*RNP1+(-1).*rcg.^(1/2).*(rcg+4.*k.*RNP1).^( ...
  1/2))+(1/2).*k.^(-1).*(rcg+2.*k.*RNP2+(-1).*rcg.^(1/2).*(rcg+4.* ...
  k.*RNP2).^(1/2))+S0)).^2)).^(1/2)).^(-1).*wc).^(1/2));
