clear, clc % 

x_bar_cg = 0.3; %
M = 484/968; % Mach
AoA_0L = -3 * (pi/180); % Radians
CMac = -0.01; %
eta_h = 0.95;


% WING %

b_w = 113; % Wingspan
cr_w = 23; % Wing Root Chord
ct_w = 4; % Wing Tip Chord
LE_sweep_w = atan(28/50); % Sweep Angle Wing in Rad
    [lambda_w, S_w, AR_w, c_bar_w, x_mgc_w, y_mgc_w] = geoParameter(b_w,cr_w,ct_w,LE_sweep_w)
    x_ac_w = c_bar_w/4;
aileronRatio = 0.2; % (cf/c)
    tau_a = AppendixD(aileronRatio);
y1 = 35; % inner point of aileron
y2 = 51; % outer point of aileron
canthalTilt = atan(9/50); % divergence - (rad)

% HORIZONTAL TAIL % 

b_h = 46; % HorTailspan
cr_h = 12; % HorTail Root Chord
ct_h = 4; % HorTail Tip Chord
LE_sweep_h = atan(13/23); % Sweep Angle HorTail in Rad
    [lambda_h, S_h, AR_h, c_bar_h, x_mgc_h, y_mgc_h] = geoParameter(b_h,cr_h,ct_h,LE_sweep_h)
elevatorRatio = 0.3;
    tau_e = AppendixD(elevatorRatio)

x_whr = 55; % dist from LE of wing root to LE of horTail chord
    x_wh = x_whr - cr_w/4 + cr_h/4;
z_wh = 6; % vert distance from wing to tail, positive if tail above
    r = x_wh/(b_w/2);
    m = z_wh/(b_w/2);
    x_ac_h = x_whr + x_mgc_h + c_bar_h/4 - x_mgc_w;
    x_cg = x_bar_cg*c_bar_w;

% VERTICAL TAIL %

b_v = 31; % dist from A/C centerline to tail tip
cr_v = 26; % root chord vertical tail
ct_v = 6; % tip chord vertical tail
eta_v = 0.95; % assumed - if not given just keep it as this
LE_sweep_v = atan(27/31); % rad
    [lambda_v, S_v, AR_v, c_bar_v, x_mgc_v, y_mgc_v] = geoParameter(b_v,cr_v,ct_v,LE_sweep_v)
    x_mgc_v=x_mgc_v*2; y_mgc_v=y_mgc_v*2;
    AR_eff = 2*AR_v;
x_wvr = 39; % distance from wing leading edge to v leading edge
    x_ac_v = x_wvr + x_mgc_v + (c_bar_v/4)-x_mgc_w;
z_w = 0; % Distance from centerline to wing % FROM THE HANDOUT LETS GO
d_max = 12; % Max diameter
rudderRatio = 0.3; 
    tau_r = AppendixD(rudderRatio)


% PARAMETER TIME :)
CL_a_w = polhamus(AR_w,lambda_w,LE_sweep_w,M);
CL_a_h = polhamus(AR_h,lambda_h,LE_sweep_h,M);
de_dalpha = downwash(AR_w, lambda_w, LE_sweep_w, M, r, m);

CL_0 = CL_a_w * abs(AoA_0L)
CL_a = CL_a_w + eta_h*(S_h)/(S_w)*(1-de_dalpha)*CL_a_h
CL_i_h = eta_h*(S_h)/(S_w)*CL_a_h
CL_d_e = CL_i_h*tau_e

Cm_0 = CMac + CL_0*(x_cg-x_ac_w)/c_bar_w
Cm_a = CL_a_w*(x_cg-x_ac_w)/c_bar_w - eta_h*(S_h/S_w)*(1-de_dalpha)*CL_a_h*(x_ac_h-x_cg)/c_bar_w
Cm_i_h = -eta_h*CL_a_h*(S_h/S_w)*(x_ac_h-x_cg)/c_bar_w
Cm_d_e = Cm_i_h*tau_e

x_bar_NP = (CL_a_w*x_ac_w + eta_h*(S_h/S_w)*(1-de_dalpha)*CL_a_h*x_ac_h)/(c_bar_w*(CL_a_w + eta_h*(S_h/S_w)*(1-de_dalpha)*CL_a_h)) 
x_NP = x_bar_NP*c_bar_w;
SM = (x_bar_NP - x_bar_cg)*100 % PERCENT

CL_a_v = polhamus(AR_eff, lambda_v, LE_sweep_v, M)
AppE = AppendixE(S_v, S_w, z_w, d_max, AR_w, LE_sweep_w, lambda_w)
Cy_0 = 0;
Cy_B_v = -S_v/S_w*CL_a_v*AppE;
Cy_B_w = -0.0001*abs(canthalTilt)*180/pi;
Cy_B_wv = Cy_B_w+Cy_B_v;
Cy_B_f = 0.3*Cy_B_wv;
Cy_B = Cy_B_f+Cy_B_wv
Cy_d_a = 0
Cy_d_r = eta_v*S_v/S_w*CL_a_v*tau_r

Cl_0 = 0;
Cl_B_w = -2*CL_a_w*canthalTilt*y_mgc_w/b_w;
Cl_B_v = Cy_B_v* y_mgc_v/b_w;
Cl_B = Cl_B_v + Cl_B_w
Cl_d_a = AppendixG(CL_a_w,tau_a,cr_w,S_w,b_w,lambda_w,y1,y2)
Cl_d_r = Cy_d_r*y_mgc_v/b_w

Cn_0 = 0;
Cn_B = -Cy_B_v*(x_ac_v-x_cg)/b_w % equal to Cn_B_v
Cn_d_a = 0
Cn_d_r = -Cy_d_r*(x_ac_v-x_cg)/b_w


function [lambda, S, AR, c_bar, x_mgc, y_mgc] = geoParameter(b,cr,ct,sweep) % various parameters that need to be found over and over
    lambda = ct/cr; % Taper Ratio
    S = 0.5*b*cr*(1+lambda); % Surface area
    AR = b^2/S; % Aspect Ratio
    c_bar = 2/3 * cr * (1+lambda+lambda^2)/(1+lambda); % self explanetory
    y_mgc = b/6 * (1+2*lambda)/(1+lambda); % same with this 
    x_mgc = y_mgc*tan(sweep); % and this as well
end

function tau = AppendixD(flapRatio) % Effectivness calculation
    tau = 1.340933 + (0.00003390316 - 1.340933)/(1 + (flapRatio/0.4437918)^1.331642); % from online curve fit, accurate except at really low ratios like <0.05
end

function CL_a = polhamus(AR,lambda,sweep,M) % CL_alpha function 
   if AR < 4
       k = 1 + (AR*(1.87-0.000233*sweep))/100;
   elseif AR >= 4
       k = 1+((8.2 -2.3*sweep)-AR*(0.22-0.153*sweep))/100;
   end
   tangent_half = tan(sweep)-(4*0.5*(1-lambda))/(AR*(1+lambda));
   CL_a = 2*pi*AR/(2+sqrt(((AR^2)*(1-M^2))/(k^2)*(1 + ((tangent_half)^2)/(1-M^2))+4));
end


function de_dalpha = downwash(AR,lambda,LE_sweep,M,r,m)
    KA = (1/AR)-(1/(1+AR^1.7));
    K_lambda = (10-3*(lambda))/7;
    Kmr = (1-(m*0.5))/(r^.33);
    Sweep_25 = atan(tan(LE_sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));
    de_o_alpha_0 = 4.44*(KA*K_lambda*Kmr*sqrt(cos(Sweep_25)))^1.19;
    de_dalpha = de_o_alpha_0/sqrt(1-M^2);
end

function AppE = AppendixE(S_v, S, z, d, AR, sweep, lambda)
    q = atan(tan(sweep)-(4*0.25*(1-lambda))/(AR*(1+lambda)));
    AppE = 0.724 + 3.06*(S_v/S)/(1+cos(q))+0.4*z/d + 0.009*AR;
end

function Cl_d_a = AppendixG(cl,t,cr,s,b,l,y1,y2)
    Cl_d_a = (2*cl*t*cr/s/b)*((y2^2 /2 + ((l-1)/(b/2))*y2^3 /3)-(y1^2 /2 + ((l-1)/(b/2))*y1^3 /3));
end
