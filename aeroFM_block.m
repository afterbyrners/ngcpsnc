function [Fx,Fy,Fz,Mx,My,Mz,alpha,beta] = aeroFM_block( ...
    u,v,w,p,q,r, delta_e,delta_a,delta_r, ...
    rho, S, b, cbar, ...
    CL0,CLa,CLde, Cm0,Cma,Cmde, ...
    Cyb, Clb, Cnb, Clda, Cydr, Cndr, ...
    CD0, OswaldE, CoeffsPerDeg)
% aeroFM_block
% this drops directly inside a Simulink MATLAB Function block.
% All inputs are scalar signals (Constant blocks for params).
%
% Outputs are scalar forces/moments (N, N*m) and angles (rad).

% Convert per-degree derivatives if flagged
if CoeffsPerDeg ~= 0
    d2r = pi/180;
    CLa  = CLa  * (180/pi);
    CLde = CLde * (180/pi);
    Cma  = Cma  * (180/pi);
    Cmde = Cmde * (180/pi);
    Cyb  = Cyb  * (180/pi);
    Clb  = Clb  * (180/pi);
    Cnb  = Cnb  * (180/pi);
    Clda = Clda * (180/pi);
    Cydr = Cydr * (180/pi);
    Cndr = Cndr * (180/pi);
end

% Kinematics
V = sqrt(u*u + v*v + w*w) + eps;
alpha = atan2(w, max(1e-9,u));
tmp = v / V;
if tmp > 0.999, tmp = 0.999; end
if tmp < -0.999, tmp = -0.999; end
beta  = asin(tmp);

% Dynamic pressure
qbar = 0.5 * rho * V*V;

% Coefficients
CL = CL0 + CLa*alpha + CLde*delta_e;
Cm = Cm0 + Cma*alpha + Cmde*delta_e;
CY = Cyb*beta + Cydr*delta_r;
Cl = Clb*beta + Clda*delta_a;
Cn = Cnb*beta + Cndr*delta_r;

% Drag model
AR = b*b / S;
if AR <= 0, AR = 1e-6; end
if OswaldE <= 0, OswaldE = 0.8; end
k  = 1.0 / (pi * OswaldE * AR);
if CD0 < 0, CD0 = 0.03; end
CD = CD0 + k * CL*CL;

% Forces
ca = cos(alpha); sa = sin(alpha);
L = qbar * S * CL;
D = qbar * S * CD;
Y = qbar * S * CY;

Fx = -(D*ca - L*sa);
Fy =  Y;
Fz = -(D*sa + L*ca);

% Moments
Mx = qbar * S * b    * Cl;
My = qbar * S * cbar * Cm;
Mz = qbar * S * b    * Cn;
end