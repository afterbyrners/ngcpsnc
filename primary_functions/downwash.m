function de_dalpha = downwash(AR,lambda,LE_sweep,M,r,m)
    KA = (1/AR)-(1/(1+AR^1.7));
    K_lambda = (10-3*(lambda))/7;
    Kmr = (1-(m*0.5))/(r^.33);
    Sweep_25 = atan(tan(LE_sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));
    de_o_alpha_0 = 4.44*(KA*K_lambda*Kmr*sqrt(cos(Sweep_25)))^1.19;
    de_dalpha = de_o_alpha_0/sqrt(1-M^2);
end
