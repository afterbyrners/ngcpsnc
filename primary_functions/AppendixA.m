% Appendix A Function

% The "Appendix A" function calculates extremely common-to-see values for a
% wing using only basic fundamental measurements.

% TODO accept inputs as taper ratio

%% INPUTS & OUTPUTS

% AR       = aspect ratio of the wing
% LE_sweep = leading edge sweep angle, in radians
% lambda   = taper ratio of the wing
% M        = flight mach number
% r        = nondimensional x-wise displacement of the surface
% m        = nondimensional z-wise displacement of the surface

% CL_a     = coefficient of lift derivative to AoA for a lifting surface

function [lambda, S, AR, c_bar, x_mgc, y_mgc] = AppendixA(b,cr,ct,sweep) % various parameters that need to be found over and over
    lambda = ct/cr; % Taper Ratio
    S = 0.5*b*cr*(1+lambda); % Surface area
    AR = b^2/S; % Aspect Ratio
    c_bar = 2/3 * cr * (1+lambda+lambda^2)/(1+lambda); % self explanetory
    y_mgc = b/6 * (1+2*lambda)/(1+lambda); % same with this 
    x_mgc = y_mgc*tan(sweep); % and this as well
end