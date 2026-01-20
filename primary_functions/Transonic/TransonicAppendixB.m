% Transonic Appendix B Function

% The "Appendix B" function calculates the variation of coefficient of lift
% with the angle of attack of the plane. 

% This range is only valid for:
% Mach 0.8 to 1.2
% Aspect Ratios between 3 and 8
% Taper Ratios between 0.4 and 10

% The transonic variation of this replaces the original Prandtl-Glauert
% term with a Karman-Tsien term.

%% INPUTS & OUTPUTS

% AR       = aspect ratio of the wing
% LE_sweep = leading edge sweep angle, in radians
% lambda   = taper ratio of the wing
% M        = flight mach number

% CL_a     = coefficient of lift derivative to AoA for a lifting surface

function CL_a = TransonicAppendixB(AR,sweep_LE,lambda,M) % CL_alpha function 
   if AR < 4
       k = 1 + (AR * (1.87 - 0.000233 * sweep_LE))/100;
   elseif AR >= 4
       k = 1 + ((8.2 - 2.3 * sweep_LE) - AR * (0.22 - 0.153 * sweep_LE)) / 100;
   end
   KT = (sqrt(1 - M^2) + (M^2 / (1 + sqrt(1 - M^2))) * 2*pi / 2);

   half_chord_sweep = atan(tan(sweep_LE)-(4*0.5*(1-lambda))/(AR*(1+lambda)));

   CL_a = 2*pi*AR / (2 + sqrt( ((AR^2)*(KT)^2)/(k^2)*(1 + (tan(half_chord_sweep)^2)/(KT)^2 )+4));
end