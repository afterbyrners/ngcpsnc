%%% FORCE BUILD-UP

% Description: This program uses the geometric and aerodynamic parameters
% of the subject aircraft in order to compute the necessary forces and
% moments needed to properly analyze the aircraft.

% Afterwards, these stability derivatives can then be implemented in a
% separate but accompanying program to assemble a perturbation apprehensive
% state-space model of the aircraft.

% It currently only works for normalized tails, meaning that it will not
% work for a plane that uses a V-tail. This is being developed further.

folder = fileparts(which("force_build_up.m"));
addpath(genpath(folder));

clear

%% FUNDAMENTAL AXIS NOTES (For a later README)

% This presumes the Stability Axis System (SAS) of the aircraft, which is
% fundamentally a modification of the Body Axis System (BAS) of the normal
% aircraft with a rotation of the X and Z axes by the trim angle of attack.

% In the BAS:
% The X-axis is +X toward the nose & -X toward the aft of the plane
% The Y-axis is +Y toward the right & -Y toward the left of the plane
% The Z-axis is +Z toward the bottom & -Z toward the top of the plane

% In SAS, the sideslip is still considered, but the perturbations in the Y
% axis remain so small that only the sideslip is measured and the velocity
% projected onto the stability axis Xs, U1s, is presumed to be the total
% airspeed of the vehicle, VP1. Also, the velocity measured on the Zs axis
% is considered to be so small, that it is not to be measured. 

%% NOTATIONS

% X  : Axis running lengthwise with the plane, +X at nose, -X at aft
% Y  : Axis running spanwise with the plane, +Y at right, -Y at left
% Z  : Axis running depthwise with the plane, +Z at bottom, -Z at top

% θ  : Pitch angle/attitude of the plane, called "theta" usually
% Φ  : Bank/roll angle/attitude of the plane, called "phi" usually
% Ψ  : Heading angle/attitude of the plane, called "psi" usually
% α  : Angle of attack of the plane, called "alpha" or "AOA"
% β  : Sideslip angle of the plane, called "beta"

% VP : True aircraft air speed 
% U  : Component of VP along X
% V  : Component of VP along Y
% W  : Component of VP along Z

% FaX: Aerodynamic force component about the X axis
% FaY: Aerodynamic force component about the Y axis
% FaZ: Aerodynamic force component about the Z axis

% l  : Rolling moment made about the X axis
% m  : Pitching moment made about the Y axis
% n  : Yawing moment made about the Z axis

% P  : Roll rate
% Q  : Pitch rate
% R  : Yaw rate

% D : Coefficient of drag
% L : Coefficient of lift
% T : Coefficient of thrust (not an aerodynamic force)

% iw : Incidence angle of wing planform
% ih : Incidence angle of horizontal flight surface
% iv : Incidence angle of vertical flight surface

% δe : Deflection of elevators 
% δa : Deflection of ailerons
% δr : Deflection of rudders
% δf : Deflection of flaps
% δv : Deflection of ruddervators
% δn : Deflection of thruster nozzles

%% Plane & Flight Characteristics

% Define Altitude
h_ft = 800 + 250; % ft
h_m = h_ft * 0.3048; % Convert altitude from feet to meters

W = 35; % weight of aircraft, lbf
g = 32.2;
m  = W / g; % mass

% Get ISA properties (SI units)
[T_K, a_mps, P_Pa, rho_kgm3] = atmosisa(h_m);

% Convert to imperial units
T_R = T_K * 9/5;               % Kelvin → Rankine
a_fps = a_mps * 3.28084;       % m/s → ft/s
P_psf = P_Pa * 0.0208854342;   % Pa → lb/ft²
rho_sl = rho_kgm3 * 0.062428;  % kg/m³ → slugs/ft³

%% Flight Conditions
M_1 = 0.0539568; % Flight Mach Number
u_1 = a_fps * M_1; % Flight forward velocity
q_1 = 0.5 * rho_sl * u_1^2; % Dynamic Pressure
eta_h = 0.95; % Horizontal Tail Dynamic Pressure Ratio
eta_v = 0.95; % Vertical Tail Dynamic Pressure Ratio

%% Airfoil Parameters (NACA 64-206)
a_L_0   = -1.0; % Zero-Lift Angle of Attack, deg
a_L_0   = deg2rad(a_L_0); % Zero-Lift Angle of Attack, rad
C_m_ac  = -0.040; % Pitching Moment at Zero-Lift Coefficient

%% Key Parameters Related to General Vehicle
X_cg = 1.63; % Actual Location of CG, ft
d = 0.90; % Diameter of fuselage, for now the largest dim, ft

% Maximum Thrust Vector Angle in both directions
dt_max_deg = 20; % deg
dt_max_rad = deg2rad(dt_max_deg); % rad

% It is worth noting that small angle approximation can be used such that
% Fz = Tsin(dt_max) = T*dt_max, since sin(dt_max) rough equals dt_max if
% dt_max is around less than 30 deg. In this sense, Fx = Tcos(dt_max) = T
% since cos(dt_max) is roughly about 0.94. % Linearizations can be made.

%% Wing Design Parameters
S_w = 8.7725; % Wing Reference Area, ft^2
AR_w = 5.58511; % Wing Aspect Ratio, N/A
sweep_w = 3.63529; % Wing Leading Edge Sweep, deg
dihedral_w = 0; % Wing Dihedral Angle, deg
taper_w = 0.55; % Wing Taper Ratio, 
cf_c_w = 0.3; % Wing Chord-Flap Ratio
x_w = 1.3; % Absolute Wing Front Point X Location
y_w = 0; % Absolute Wing Front Point Y Position
z_w = 0.315; % Absolute Wing Front Point Z Position

% Convert Wing Leading Edge Sweep to radians
sweep_w = deg2rad(sweep_w); % rad

% Calculate Wing Span
b_w = sqrt(S_w * AR_w); % ft

% Calculate Wing Root Chord
c_r_w = (2 * S_w) / (b_w * (1 + taper_w)); % ft

% Calculate Wing Tip Chord
c_t_w = c_r_w * taper_w; % ft

% Calculate Wing Mean Geometric Chord (MGC)
mgc_w = 2/3 * c_r_w * (1 + taper_w + taper_w ^ 2) / (1 + taper_w); % ft

% Calculate Y-position of Wing MGC
y_mgc_w = (b_w / 6) * (1 + 2 * taper_w) / (1 + taper_w); % ft

% Calculate X-position of Wing MGC
x_mgc_w = y_mgc_w * tan(sweep_w); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_w = RaymerWingLift(AR_w,sweep_w,taper_w,M_1); % 1/rad

% Calculate Effectiveness of Wing Control Surfaces (Ailerons)
tau_w = AppendixD(cf_c_w); % N/A

% Calculate Non-Dimensional Wing Aerodynamic Center Position X-Wise
x_ac_w = ((x_w + x_mgc_w + 0.25 * mgc_w) - (x_w + x_mgc_w)) / mgc_w; % N/A

% Calculate Non-Dimensional Center of Gravity Position X-Wise
x_cg = ((X_cg) - (x_w + x_mgc_w)) / mgc_w; % N/A

%% Horizontal Tail (H.T.) Design Parameters
S_h = 1.69; % H.T. Reference Area, ft^2
AR_h = 3; % H.T. Aspect Ratio, N/A
sweep_h = 20; % H.T. Leading Edge Sweep, deg
taper_h = 0.5; % H.T. Taper Ratio, 
cf_c_h = 0.3; % H.T. Chord-Flap Ratio
x_h = 4.9; % Absolute H.T. Front Point X Location
y_h = 0; % Absolute H.T. Front Point Y Position
z_h = 0.2; % Absolute H.T. Front Point Z Position

% Convert leading edge sweep to radians
sweep_h = deg2rad(sweep_h); % rad

% Calculate Wing Span
b_h = sqrt(S_h * AR_h); % ft

% Calculate Root Chord
c_r_h = (2 * S_h) / (b_h * (1 + taper_h)); % ft

% Calculate Tip Chord
c_t_h = c_r_h * taper_h; % ft

% Calculate Mean Geometric Chord (MGC)
mgc_h = 2/3 * c_r_h * (1 + taper_h + taper_h ^ 2) / ... 
    (1 + taper_h); % ft

% Calculate Y-position of MGC
y_mgc_h = (b_h / 6) * (1 + 2 * taper_h) / (1 + taper_h); % ft

% Calculate X-position of MGC
x_mgc_h = y_mgc_h * tan(sweep_h); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_h = RaymerWingLift(AR_h,sweep_h,taper_h,M_1); % 1/rad

% Calculate Effectiveness of H.T. Control Surfaces (Elevator)
tau_h = AppendixD(cf_c_h); % N/A

% Calculate Non-Dimensional H.T. Aerodynamic Center Position X-Wise
x_ac_h = ((x_h + x_mgc_h + 0.25 * mgc_h) - (x_w + x_mgc_w)) / mgc_w; % N/A

% Calculate X Distance Between ACs of H.T. and Wing
x_w_h = x_h - (c_r_w / 4) + (c_r_h / 4);

% Calculate Non-Dimensional X Distance between H.T. and Wing
r_h = x_w_h / (0.5 * b_w);

% Calculate Z Distance Between ACs of Horizontal Tail and Wing 
z_w_h = (z_h) - z_w;

% Calculate Non-Dimensional Z Distance between H.T. and Wing
m_h = z_w_h / (0.5 * b_w);

% Calculate downwash gradient for horizontal tail equivalent
dw_h = AppendixC(AR_w,sweep_w,taper_w,M_1,r_h,m_h);

%% Vertical Tail (V.T.) Design Parameters
S_v = 0.78610; % V.T. Reference Area, ft^2
AR_v = 2; % V.T. Aspect Ratio, N/A
sweep_v = 20; % V.T. Leading Edge Sweep, deg
taper_v = 0.5; % V.T. Taper Ratio, 
cf_c_v = 0.3; % V.T. Chord-Flap Ratio
x_v = 5; % Absolute V.T. Front Point X Location
y_v = 0; % Absolute V.T. Front Point Y Position
z_v = 0.167; % Absolute V.T. Front Point Z Position

% Convert leading edge sweep to radians
sweep_v = deg2rad(sweep_v); % rad

% Calculate Wing Span
b_v = sqrt(S_v * AR_v); % ft

% Calculate Root Chord
c_r_v = (2 * S_v) / (b_v * (1 + taper_v)); % ft

% Calculate Tip Chord
c_t_v = c_r_v * taper_v; % ft

% Calculate Mean Geometric Chord (MGC)
mgc_v = 2/3 * c_r_v * (1 + taper_v + taper_v ^ 2) / ... 
    (1 + taper_v); % ft

% Calculate Y-position of MGC
z_mgc_v = (b_v / 6) * (1 + 2 * taper_v) / (1 + taper_v); % ft

% Calculate X-position of MGC
x_mgc_v = z_mgc_v * tan(sweep_v); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_v = RaymerWingLift(AR_v,sweep_v,taper_v,M_1); % 1/rad

% Calculate Effectiveness of V.T. Control Surfaces (Elevator)
tau_v = AppendixD(cf_c_v); % N/A

% Calculate Non-Dimensional V.T. Aerodynamic Center Position X-Wise
x_ac_v = ((x_v + x_mgc_v + 0.25 * mgc_v) - (x_w + x_mgc_w)) / mgc_w; % N/A

% Calculate X Distance Between ACs of V.T. and Wing
x_w_v = x_v - (c_r_w / 4) + (c_r_v / 4);

% Calculate Z Distance Between ACs of V.T. and Wing 
z_w_v = (z_v + z_mgc_v) - z_w;

%% D Forces (Drag)

% Calculate Trim Drag Coefficient
C_D_1 = 0.03; % This value makes sense for an f-22. Usually CDx1 but changed to CD1 for experimental

% Calculate Variation of Drag Coefficient with Angle of Attack
C_D_a = 0.07;

% Variation of airplane drag coefficient with dimensionless speed
C_D_u = 0; % Assumed

% Variation of airplane drag coefficient with elevator deflection angle
C_D_de = 0; % Assumed

%% T Forces (Thrust)

% Calculate Trim Thrust Coefficient
C_T_x_1 = C_D_1;

% Calculate Variation of Thrust Coefficient with Fwd V Perturbation
C_T_x_u = M_1 / (q_1 * S_w) - (2 * C_T_x_1); % Roskam Equation 3.206

% Calculate Variation of X Thrust with Nozzle Deflection
C_T_x_d_n = C_T_x_1 * (cos(dt_max_rad)) / dt_max_rad;

% Calculate Variation of Z Thrust with Nozzle Deflection
C_T_z_d_n = C_T_x_1 * (sin(dt_max_rad)) / dt_max_rad;

%% L Forces (Lift)

% Calculate Level Flight Coefficient of Lift
%C_L_1 = -(C_L_a_w * a_L_0);  % N/A
C_L_1 = W / (q_1 * S_w); %Double check this with zero-lift airfoil versus trim

% Calculate Coefficient of Lift caused by General Plane
C_L_a = C_L_a_w + eta_h * (S_h / S_w) * C_L_a_h * (1 - dw_h);  % 1/rad

% Calculate Coefficient of Lift caused by Horizontal Stabilizer Incidence
C_L_ih = eta_h * (S_h / S_w) * C_L_a_h;  % 1/rad

% Calculate Coefficient of Lift caused by Elevator Flaps
C_L_de = C_L_ih * tau_h;  % 1/rad

% Calculate Coefficient of Lift caused by pitch rate
C_L_q = 2 * C_L_a_h * (x_ac_h - x_cg) * (S_h/S_w);

% Calculate Coefficient of Lift caused by change in dimensionless x speed u
C_L_u = 0; % Assumed

%% X Forces (Drag & Thrust)

% Will change later

%% Calculate Y Forces (Sideforces)

% Calculate Level Flight Coefficient of Side Force
C_Y_0 = 0;  % Given/assumed due to symmetry
 
% Calculate Coefficient of Side Force caused by V.T. sideslip
C_Y_b_v = -(S_v / S_w) * C_L_a_v * AppendixE(S_v, S_w, z_w, d, AR_w, sweep_w, taper_w);

% Calculate Coefficient of Side Force caused by wing sideslip 
C_Y_b_w = -0.0001 * abs(dihedral_w) * 180 / pi;  % rad^-1, basically zero

% Calculate Coefficient of Side Force caused by fuselage sideslip
C_Y_b_f = 0.3 * C_Y_b_w;

% Calculate Coefficient of Side Force caused by sideslip
C_Y_b = C_Y_b_f + C_Y_b_v + C_Y_b_w;

% Calculate change in Coefficient of Side Force caused by aileron incidence
C_Y_da = 0;  % we can assume so due to ailerons not affecting that DOF

% Calculate change in Coefficient of Side Force caused by rudder incidence
C_Y_dr = eta_v * (S_v / S_w) * C_L_a_v * tau_v;

% Calculate change in Coefficient of Side Force caused by roll rate
C_Y_p = -2 * C_L_a_v * (z_w_v / b_w) * eta_v * (S_v/S_w);

% Calculate change in Coefficient of Side Force caused by yaw rate
C_Y_r = C_L_a_v * (2 * ((x_ac_v - x_cg) * mgc_w)) / b_w * eta_v * (S_v/S_w);

%% Calculate Z Forces (Lift(?))

% TODO

%% Calculate l Moments (Rolling)

% Calculate Level Flight Coefficient of Rolling Moment
C_l_0 = 0;  % Assumed due to symmetry

% Calculate Coefficient of Rolling Moment caused by V.T. sideslip
C_l_b_v = C_Y_b_v * (z_w_v / b_w);

% Calculate Coefficient of Rolling Moment caused by wing sideslip
C_l_b_w = -2 * C_L_a_w * dihedral_w * (y_mgc_w / b_w);

% Calculate Coefficient of Rolling moment caused by sideslip
C_l_b  = C_l_b_v + C_l_b_w;

% Calculate Coefficient of Rolling Moment change caused by ailerons
a_1 = 0.5 * (0.5 * b_w);
a_2 = 1.0 * (0.5 * b_w);

C_l_da = AppendixG(C_L_a_w,tau_w,c_r_w,S_w,b_w,taper_w,a_1,a_2);

% Calculate Coefficient of Rolling Moment change caused by rudder
C_l_dr = C_Y_dr * (z_w_v / b_w);

% Calculate Coefficient of Rolling moment caused by yaw rate 
C_l_r = C_L_a_v * 2 * ((x_ac_v - x_cg) * mgc_w) * z_w_v / b_w^2 * eta_v ...
        * (S_v/S_w); 
% Consider finding the fuselage contributions.

% Calculate Coefficient of Rolling Moment caused by roll rate of V.T.
C_l_p_v = -2 * C_L_a_v * (z_w_v / b_w)^2 * eta_v * (S_v/S_w); 
% ^ Consider finding the fuselage & htail contributions.

% Calcualte Coefficient of Rolling Moment caused by roll rate of wing
C_l_p_wf = -(C_L_a_w / 4) * AR_w;

% Calculate Coefficient of Rolling Moment caused by roll rate
C_l_p = C_l_p_wf + C_l_p_v;

%% Calculate m Moments (Pitching)

% Calculate Level Flight Coefficient of Moment
C_m_0 = C_m_ac + C_L_1 * (x_cg - x_ac_w);  % N/A

% Calculate Coefficient of Moment caused by General Plane
C_m_a = C_L_a_w * (x_cg - x_ac_w) - C_L_a_h * eta_h * (S_h / S_w) * ...
        (1 - dw_h) * (x_ac_h - x_cg); % 1/rad

% Calculate Coefficient of Moment caused by Horizontal Stabilizer Incidence
C_m_ih = -C_L_a_h * eta_h * (S_h / S_w) * (x_ac_h - x_cg); % 1/rad

% Calculate Coefficient of Moment caused by Flap Angles
C_m_de = -C_L_a_h * eta_h * (S_h / S_w) * tau_h * (x_ac_h - x_cg); % 1/rad

% Calculate Coefficient of Moment caused by pitch rate
C_m_q = -2.2 * C_L_a_h * eta_h * (x_ac_h - x_cg)^2 * (S_h/S_w);

% Variation of pitching moment coefficient with dimensionless speed
C_m_u = 0; % Assumed, change later

% Steady state pitching moment coefficient
C_m_1 = 0; % Assumed, change later

%% Calculate n Moments (Yawing)

% Calculate Level Flight Coefficient of Yawing Moment
C_n_0 = 0;  % assumed

% Calculate Coefficient of Yawing Moment change caused by beta
C_n_b = -C_Y_b_v * ((x_ac_v - x_cg) * mgc_w) / b_w;

% Calculate Coefficient of Yawing Moment change caused by ailerons
C_n_da = 0;  % assumed so (low dihedral = min effect)

% Calculate Coefficient of Yawing Moment change caused by rudder
C_n_dr = -C_Y_dr * ((x_ac_v - x_cg) * mgc_w) / b_w;

% Calculate Coefficient of Yawing Moment change due to thrust from sideslip
C_n_T_b = 0; % Assumed because of symmetry

% Calculate Coefficient of Yawing Moment change caused by roll rate
C_n_p = 2 * C_L_a_v * z_w_v * ((x_ac_v - x_cg) * mgc_w) * z_w_v / b_w^2 ...
        * eta_v * (S_v/S_w);
% ^ It is mentioned to be negligible, but consider fuselage contributions.

% Calculate Coefficient of Yawing Moment change caused by yaw rate
C_n_r = -C_L_a_v * 2 * ((x_ac_v - x_cg) * mgc_w)^2 / b_w^2 * eta_v * ...
        (S_v/S_w);
% ^ As with others, this contribution dominates but consider finding wing
% and fuselage contributions.

%% Calculate Flying Qualities

% Calculate Neutral Point
x_np = (x_ac_w + ((C_L_a_h / C_L_a_w) * eta_h * (S_h / S_w) * ... 
       (1 - dw_h) * x_ac_h)) / (1 + (C_L_a_h / C_L_a_w) * ... 
       eta_h * (S_h / S_w) * (1 - dw_h));

X_np = x_np * mgc_w + x_mgc_w + x_w;

% Calculate Static Margin Percentage
SM = (x_np - x_cg) * 100;

%% Important Inertia Values
Ixx = 5;
Iyy = 10;
Izz = 12;

%% State Space Coefficients: X-Axis Components (Forward-Backward)

% Change in X position due to changes in forward airspeed
Xu  = -(q_1 * S_w) / (m * u_1) * (C_D_u + 2 * C_D_1);

% Change in X position due to changes in thrust caused by speed changes
XTu = (q_1 * S_w) / (m * u_1) * (C_T_x_u + 2 * C_T_x_1);

% Change in X position due to changes in angle of attack
Xa  = (q_1 * S_w) / (m) * (-C_D_a + C_L_1);

% Change in X position due to changes in elevator deflection
Xde = (q_1 * S_w) / (m) * (C_D_de);

%% State Space Coefficients: Z-Axis Components (Downward-Upward)

% Change in Z position due to changes in forward airspeed
Zu  = -(q_1 * S_w) / (m * u_1) * (C_L_u + 2 * C_L_1);

% Change in Z position due to changes in angle of attack
Za  = -(q_1 * S_w) / (m) * (C_L_a + C_D_1);

% Variation of Z position with pitch rate
Zq  = -(q_1 * S_w * mgc_w) / (2 * m * u_1) * (C_L_q);

% Variation of Z position with elevator deflection
Z_de = -(q_1 * S_w) / (m) * (C_L_de);

%% State Space Coefficients: M Components (Pitching Moment About Y)

% Variation of pitching moment due to changes in forward airspeed
Mu  = (q_1 * S_w * mgc_w) / (Iyy * u_1) * (C_m_u + 2 * C_m_1);

% Variation of pitching moment due to changes in angle of attack
Ma  = (q_1 * S_w * mgc_w) / (Iyy) * (C_m_a);

% Variation of pitching moment due to changes in pitch rate
Mq  = (q_1 * S_w * mgc_w^2) / (2 * Iyy * u_1) * (C_m_q);

% Variation of pitching moment due to changes in pitch rate
Mde  = (q_1 * S_w * mgc_w) / (Iyy) * (C_m_de);

%% State Space Coefficients: Y Components (Right-Left)

% Variation of Y position with sideslip angle
Yb = (q_1 * S_w) / (m) * (C_Y_b);

% Variation of Y position with roll rate
Yp = (q_1 * S_w) / (2 * m * u_1) * (C_Y_p);

% Variation of Y position with rudder deflection
Ydr = (q_1 * S_w) / (m) * (C_Y_dr);

%% State Space Coefficients: L Components (Rolling Moment About X)

% Variation of rolling moment with sideslip angle
Lb  = (q_1 * S_w * b_w) / (Ixx) * (C_l_b);

% Variation of rolling moment with roll rate
Lp  = (q_1 * S_w * b_w^2) / (2 * Ixx * u_1) * (C_l_p);

% Variation of rolling moment with yaw rate
Lr  = (q_1 * S_w * b_w^2) / (2 * Ixx * u_1) * (C_l_r);

% Variation of rolling moment with aileron deflection
Lda = (q_1 * S_w * b_w) / (Ixx) * (C_l_da);

% Variation of rolling moment with rudder deflection
Ldr = (q_1 * S_w * b_w) / (Ixx) * (C_l_dr);

%% State Space Coefficients: N Components (Yawing Moment About Z)

% Variation of yawing moment with sideslip angle
Nb  = (q_1 * S_w * b_w) / (Izz) * (C_n_b);

% Variation of yawing moment with roll rate
Np  = (q_1 * S_w * b_w^2) / (2 * Izz * u_1) * (C_n_p);

% Variation of yawing moment with yaw rate
Nr  = (q_1 * S_w * b_w^2) / (2 * Izz * u_1) * (C_n_r);

% Variation of yawing moment with aileron deflection
Nda = (q_1 * S_w * b_w) / (Izz) * (C_n_da);

% Variation of yawing moment with rudder deflection
Ndr = (q_1 * S_w * b_w) / (Izz) * (C_n_dr);

%% Longitudinal State Space Assembly

% Recursive Variables:
% airspeed (U)
% angle of attack (AOA)
% pitch rate (q)
% and pitch angle (theta)

% Control Input(s): 
% elevator deflection (delta e)

% Plant Matrix Shell
A_Long      = zeros(4,4);

% First Column Inputs (U)
A_Long(1,1) = Xu + XTu; 
A_Long(2,1) = Zu / u_1; 
A_Long(3,1) = Mu;
A_Long(4,1) = 0;

% Second Column Inputs (AOA)
A_Long(1,2) = Xa;
A_Long(2,2) = Za / u_1;
A_Long(3,2) = Ma;
A_Long(4,2) = 0;

% Third Column Inputs (q)
A_Long(1,3) = 0;
A_Long(2,3) = 1 + (Zq / u_1);
A_Long(3,3) = Mq;
A_Long(4,3) = 1;

% Fourth Column Inputs (theta)
A_Long(1,4) = -g;
A_Long(2,4) = 0;
A_Long(3,4) = 0;
A_Long(4,4) = 0;

% Control Matrix Shell
B_Long      = zeros(4,1);

% First Column Inputs (Elevator deflection)
B_Long(1,1)   = Xde; 
B_Long(2,1)   = Z_de / u_1; 
B_Long(3,1)   = Mde; 
B_Long(4,1)   = 0;

%% Lateral State Space Assembly

% Recursive Variables:
% Change in sideslip (beta)
% Change in roll rate (p)
% Change in yaw rate (r)
% Change in roll angle (phi)

% Control Input(s): 
% Aileron deflection
% Rudder deflection

% Plant Matrix Shell
A_Late      = zeros(4,4);

% First Column Inputs (Beta)
A_Late(1,1) = Yb/u_1;
A_Late(2,1) = Lb; 
A_Late(3,1) = Nb;
A_Late(4,1) = 0;

% Second Column Inputs (p)
A_Late(1,2) = Yp/u_1;
A_Late(2,2) = Lp;
A_Late(3,2) = Np;
A_Late(4,2) = 1;

% Third Column Inputs (r)
A_Late(1,3) = -1;
A_Late(2,3) = Lr;
A_Late(3,3) = Nr;
A_Late(4,3) = 0;

% Third Column Inputs (r)
A_Late(1,4) = g/u_1;
A_Late(2,4) = 0;
A_Late(3,4) = 0;
A_Late(4,4) = 0;

% Control Matrix Shell
B_Late      = zeros(4,2);

% First Column Inputs (Aileron Deflection)
B_Late(1,1) = 0; 
B_Late(2,1) = Lda; 
B_Late(3,1) = Nda; 
B_Late(4,1) = 0;

% Second Column Inputs (Rudder Deflection)
B_Late(1,2) = Ydr / u_1; 
B_Late(2,2) = Ldr; 
B_Late(3,2) = Ndr; 
B_Late(4,2) = 0;

%% Directional Dynamics Analysis

% Longitudinal Dynamics
[V_Long,G_Long] = eig(A_Long);
eig_Long = diag(G_Long); % Gives Poles
tau_Long = -1 * 1 ./ (real(eig_Long)); % Time constants
w_n_Long = abs(eig_Long); % Natural frequencies
zeta_Long = -1 * real(eig_Long) ./ w_n_Long; % Damping ratios

% Short Period Calculations
[w_n_sh, idx_sh]  = max(w_n_Long);
zeta_sh = zeta_Long(idx_sh);
tau_sh  = tau_Long(idx_sh);

% Phugoid Motion Calculations
[w_n_ph, idx_ph]  = min(w_n_Long);
zeta_ph = zeta_Long(idx_ph);
tau_ph  = tau_Long(idx_ph);

% Lateral Dynamics
[V_Late,G_Late] = eig(A_Late);
eig_Late = diag(G_Late); % Gives Poles

tau_Late = -1 ./ (real(eig_Late)); % Time constants
w_n_Late = abs(eig_Late); % Natural frequencies
zeta_Late = -1 * real(eig_Late) ./ w_n_Late; % Damping ratios

% Spiral Motion Calculations
[tau_sp, idx_sp]  = max(tau_Late);
zeta_sp = zeta_Late(idx_sp);
w_n_sp  = w_n_Late(idx_sp);

% Rolling Motion Calculations
[tau_ro, idx_ro]  = min(tau_Late);
zeta_ro = zeta_Late(idx_ro);
w_n_ro  = w_n_Late(idx_ro);

% Dutch Roll Calculations
[zeta_dr, idx_dr]  = min(zeta_Late);
tau_dr = tau_Late(idx_dr);
w_n_dr  = w_n_Late(idx_dr);
