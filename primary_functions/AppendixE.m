% Appendix E Function

% The "Appendix E" function calculates the sidewash gradient of the subject
% vertical tail, based on the characteristics of the wing that it is
% dealing with and the dimensions of the fuselage of the plane.

%% INPUTS & OUTPUTS

% S_vert  = area of the vertical tail
% S_wing  = area of the reference wing
% Z_wing  = depthwise displacement of the wing
% d       = maximum fuselage depth
% AR_wing = aspect ratio of the reference wing

% AppendixE = vertical tail effectiveness, multiplied by one plus the
% sidewash gradient power.

function AppendixE = AppendixE(S_v, S, z, d, AR, sweep, lambda)
    q = atan(tan(sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));
    AppendixE = 0.724 + 3.06*(S_v/S)/(1+cos(q))+0.4*z/d + 0.009*AR;
end