% Appendix G Function

% The "Appendix G" function calculates the rolling moment contributions 
% caused by the opposite but symmetric deflections of the ailerons on the 
% main wing of the aircraft. 

%% INPUTS & OUTPUTS

% cl   = derivative of wing lift in regards to angle of attack
% t    = effectiveness multiplier of the ailerons
% cr   = root chord of the wing
% S    = area of the wing
% y1   = the inboard station spanwise loc of the aileron (closer to root)
% y2   = the outboard station spanwise loc of the aileron (closer to tip)
% l    = taper ratio of the wing   

% Cl_d_a = derivative of rolling moment caused by the aileron versus AoA

function Cl_d_a = AppendixG(cl,t,cr,s,b,l,y1,y2)
    Cl_d_a = (2*cl*t*cr/s/b)*((y2^2 /2 + ((l-1)/(b/2))*y2^3 /3)-(y1^2 /2 + ((l-1)/(b/2))*y1^3 /3));
end
