function [b_w,c_r_w,c_t_w,mgc_w,y_mgc_w,x_mgc_w,C_L_a_w,tau_w] = primaryWing(S_w,AR_w,sweep_w,taper_w,cf_c_w,M)
% Primary Wing
% This function calculates the necessary Roskam parameters for the main
% wing of the plane, presuming the plane isn't a tandem wing.

% Convert leading edge sweep to radians
sweep_w = deg2rad(sweep_w); % rad

% Calculate Wing Span
b_w = sqrt(S_w * AR_w); % ft

% Calculate Root Chord
c_r_w = (2 * S_w) / (b_w * (1 + taper_w)); % ft

% Calculate Tip Chord
c_t_w = c_r_w * taper_w; % ft

% Calculate Mean Geometric Chord (MGC)
mgc_w = 2/3 * c_r_w * (1 + taper_w + taper_w ^ 2) / ... 
    (1 + taper_w); % ft

% Calculate Y-position of MGC
y_mgc_w = (b_w / 6) * (1 + 2 * taper_w) / (1 + taper_w); % ft

% Calculate X-position of MGC
x_mgc_w = y_mgc_w * tan(sweep_w); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_w = AppendixB(AR_w,sweep_w,taper_w,M); % 1/rad

% Calculate Effectiveness of Wing Control Surfaces (Ailerons)
tau_w = AppendixD(cf_c_w); % N/A

end