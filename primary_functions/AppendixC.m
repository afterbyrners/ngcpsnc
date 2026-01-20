% Appendix C Function

% The "Appendix C" function calculates the wash gradient based on both the
% experienced Mach number of the aircraft, and the characteristics of the
% reference wing. It accounts for both up and down wash, and it is mostly
% used for the horizontal tail, unless there is an exotic tail
% configuration being used such as a V-tail or a canard.

%% INPUTS & OUTPUTS

% AR       = aspect ratio of the wing
% LE_sweep = leading edge sweep angle, in rad
% lambda   = taper ratio of the wing
% M        = flight mach number
% r        = nondimensional x-wise displacement of the surface
% m        = nondimensional z-wise displacement of the surface

% de_da    = wash gradient


function de_dalpha = AppendixC(AR,LE_sweep,lambda,M,r,m)

    KA = (1/AR)-(1/(1+AR^1.7));

    K_lambda = (10-3*(lambda))/7;

    Kmr = (1-abs(m*0.5))/(r^.33);

    Sweep_25 = atan(tan(LE_sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));

    de_o_alpha_0 = 4.44*(KA*K_lambda*Kmr*sqrt(cos(Sweep_25)))^1.19;

    de_dalpha = de_o_alpha_0/sqrt(1-M^2);

end

% Investigate the absolute value edit