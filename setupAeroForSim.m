% setupAeroForSim.m
% Build aero model from  wrapped function, convert geometry to SI,
% and publish constants to the base workspace for Simulink.
%% 1) flight-condition / model params
p = struct();
p.x_bar_cg   = 0.30;      % CG as fraction of MAC
p.M          = 484/968;   % Mach
p.AoA_0L_deg = -3;        % alpha for CL = 0 (deg)
p.CMac       = -0.01;     % Cm about MAC at CL=0
p.eta_h      = 0.95;      % tail efficiency

%% 2) build aero model  (original code wrapped as a function)
aero = buildAeroModel(p);

%% 3) UNIT CONSISTENCY: convert geometry to SI (meters, m^2)
% geometry coming out of buildAeroModel is in inches/in^2 (based on magnitudes).
% Convert ONLY the reference geometry that gets used to dimensionalize forces/moments.
inch_to_m = 0.0254;
in2_to_m2 = inch_to_m^2;

aero.S    = aero.S    * in2_to_m2;   % in^2 -> m^2
aero.b    = aero.b    * inch_to_m;   % in   -> m
aero.cbar = aero.cbar * inch_to_m;   % in   -> m

%% 4) publish geometry to base workspace (Simulink Constant blocks read these)
assignin('base','S',    aero.S);
assignin('base','b',    aero.b);
assignin('base','cbar', aero.cbar);

%% 5) publish longitudinal derivatives
assignin('base','CL0',  aero.CL0);
assignin('base','CLa',  aero.CLa);
assignin('base','CLde', aero.CLde);
assignin('base','Cm0',  aero.Cm0);
assignin('base','Cma',  aero.Cma);
assignin('base','Cmde', aero.Cmde);

%% 6) publish lateral-directional derivatives
assignin('base','Cyb',  aero.Cyb);
assignin('base','Clb',  aero.Clb);
assignin('base','Cnb',  aero.Cnb);
assignin('base','Clda', aero.Clda);
assignin('base','Cydr', aero.Cydr);
assignin('base','Cndr', aero.Cndr);

%% 7) environment & mass properties (SI)
assignin('base','rho', 1.225);  % kg/m^3
assignin('base','m',   8.0);    % kg   (replace with actual mass)
assignin('base','Ix',  0.35);   % kg·m^2  (replace with inertias)
assignin('base','Iy',  0.50);
assignin('base','Iz',  0.65);

%% 8) (optional) quick echo so you can see what's going in
disp('--- setupAeroForSim: key values (SI) ---');
fprintf('S = %.4f m^2,  b = %.4f m,  cbar = %.4f m\n', aero.S, aero.b, aero.cbar);
fprintf('CL0=%.4f  CLa=%.4f  CLde=%.4f\n', aero.CL0, aero.CLa, aero.CLde);
fprintf('Cm0=%.4f  Cma=%.4f  Cmde=%.4f\n', aero.Cm0, aero.Cma, aero.Cmde);
fprintf('Cyb=%.4f  Clb=%.4f  Cnb=%.4f  Clda=%.4f  Cydr=%.4f  Cndr=%.4f\n', ...
    aero.Cyb, aero.Clb, aero.Cnb, aero.Clda, aero.Cydr, aero.Cndr);
disp('-----------------------------------------');
