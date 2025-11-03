function AppE = AppendixE(S_v, S, z, d, AR, sweep, lambda)
    q = atan(tan(sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));
    AppE = 0.724 + 3.06*(S_v/S)/(1+cos(q))+0.4*z/d + 0.009*AR;
end