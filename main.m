%% Constants
clear, clc
v = 456; % ft/sec
h = 30000; % ft
M = 0.459; % mach
W = 6360; % lbs
q_bar = 92.7; % dynamic pressure
c_bar = 5.47; % c_bar
AoA = deg2rad(2); % angle of attack
IxxB = 7985; % slug*ft^2
IyyB = 3326; % slug*ft^2
IzzA = 11183; % slug*ft^2
IxzB = 0; % slug*ft^2
g = 32.2;
S = 182; 
m = W/g;
b = 33.8;
U1 = v;
%% Given Derivatives & Calculated Derivatives
Ctxu=-0.07;
Ctx1=0.0300;
Cdu=0;
Cd1=0.0300;
Cdalpha=0.250;
Cddeltae=0;
Clu=0;
Clalpha=5.15;
Cl1=0.378;
Cldeltae=0.5;
Clq=4.1;
Cmu=0;
Cm1=0;
Cmalpha=-0.700;
Cmq=-14.9;
Cmdeltae=-1.12;
Xu=-q_bar*S/(m*U1)*(Cdu+2*Cd1);
XTu=q_bar*S/(m*U1)*(Ctxu+2*Ctx1);
Xalpha=q_bar*S/(m)*(-Cdalpha+Cl1);
Xdeltae=-q_bar*S/m*Cddeltae;
Zalpha=-q_bar*S/m*(Clalpha+Cd1);
Zw = Zalpha/U1;
Zu=-q_bar*S/(m*U1)*(Clu+2*Cl1);
Zq=-q_bar*S*c_bar*Clq/(2*m*U1);
Zdeltae=-q_bar*S/m*Cldeltae;
Mu=q_bar*S*c_bar*(Cmu+2*Cm1)/(IyyB*U1);
Malpha=q_bar*S*c_bar*Cmalpha/IyyB;
Mq=q_bar*S*c_bar^2*Cmq/(2*IyyB*U1);
Mdeltae=q_bar*S*c_bar*Cmdeltae/IyyB;
CyB = -.346;
Cyp = -.0827;
Cyr = 0.3;
Cydeltaa = 0;
Cydeltar = 0.2;
ClB = -0.0944;
Clp = -0.442;
Clr = .0926;
Cldeltaa = 0.181;
Cldeltar = 0.015;
CnB = .1106;
Cnp = -.0243;
Cnr = -.139;
Cndeltaa = -.0254;
Cndeltar = -.0365;
YB = q_bar*S/m*CyB;
YP = q_bar*S/(2*m*U1)*Cyp;
Yr = q_bar*S/(2*m*U1)*Cyr;
Ydeltaa = 0;
Ydeltar = q_bar*S/m*Cydeltar;
LB = q_bar*S*b/IxxB*ClB;
LP = q_bar*S*b^2/(2*IxxB*U1)*Clp;
Lr = q_bar*S*b^2/(2*IxxB*U1)*Clr;
Ldeltaa = q_bar*S*b/IxxB*Cldeltaa;
Ldeltar = q_bar*S*b/IxxB*Cldeltar;
NB = q_bar*S*b/IzzA*CnB;
NP = q_bar*S*b^2/(2*IzzA*U1)*Cnp;
Nr = q_bar*S*b^2/(2*IzzA*U1)*Cnr;
Ndeltaa = q_bar*S*b/IzzA*Cndeltaa;
Ndeltar = q_bar*S*b/IzzA*Cndeltar;
%% Longitudinal
A = [(Xu+XTu) Xalpha 0 -g;
Zu/U1 Zw Zq/U1+1 0;
Mu Malpha Mq 0;
0 0 1 0];
B = [Xdeltae Zdeltae/U1 Mdeltae 0]';
C = eye(4);
D = zeros(4,1);
sys_lon = ss(A,B,C,D);

% MODES
wn_sp = sqrt(Zw*Mq-Malpha); %fprintf("\nThe Short Period Mode's Natural Frequency is %.2f\n",wn_sp)
zeta_sp = -(Mq +Zw)/(2*wn_sp); %fprintf("\nThe Short Period Mode's Damping Ratio is %.2f\n",zeta_sp)
tau_sp = 1/(wn_sp*zeta_sp); %fprintf("\nThe Short Period's Time Constant is %.2f\n",tau_sp)
wn_ph = sqrt(-g*Zu/U1); %fprintf("\nThe Phugoid Mode's Natural Frequency is %.2f\n",wn_ph)
zeta_ph = -(Xu+XTu)/(2*wn_ph); %fprintf("\nThe Phugoid Mode's Damping Ratio is %.2f\n",zeta_ph)
tau_ph = 1/(wn_ph*zeta_ph); %fprintf("\nThe Phugoid's Time Constant is %.2f\n",tau_ph)

% real values for modes:
poles_lon = eig(A); % finds poles (short period is first two and ph is bottom two)
zeta_lon_real = abs(cos(angle(poles_lon))); % matlab website is so goated ^o^
wn_lon_real = abs(real(poles_lon))./zeta_lon_real; % 1 and 2 = sp, 3 and 4 = ph
tau_lon_real = 1./(abs(real(poles_lon))); % same order as other ones

% problem 5
[b_lon,a_lon] = ss2tf(A,B,C,D,1);
% DU_DE_tf = tf(U1*b_lon(1,:),U1*a_lon)
DA_DE_tf = tf(U1*b_lon(2,:),U1*a_lon); % AoA per elevator
% DQ_DE_tf = tf(U1*b_lon(3,:),U1*a_lon) 
DTHE_DE_tf = tf(U1*b_lon(4,:),U1*a_lon); % pitch angle per elevator

figure, rlocus(DA_DE_tf), title("AoA to Elevator Transfer Function");
figure, rlocus(DTHE_DE_tf), title("Pitch Angle to Elevator Transfer Function");
figure, bode(DA_DE_tf), margin(DA_DE_tf), title("AoA to Elevator Transfer Function")
figure, bode(DTHE_DE_tf), margin(DTHE_DE_tf), title("Pitch Angle to Elevator Transfer Function");
%% Lateral
A = [YB/U1, YP/U1, -1, g/U1;
LB, LP, Lr, 0;
NB, NP, Nr, 0;
0, 1, 0, 0];
B = [0 Ldeltaa Ndeltaa 0;
Ydeltar/U1 Ldeltar Ndeltar 0]';
C = eye(4);
D = zeros(4,2);
sys_lat = ss(A,B,C,D);

% MODES
zeta_r = 1; %fprintf("\nThe Roll Mode's Damping Ratio is %.2f\n",zeta_r)
tau_r = -1/LP; %fprintf("\nThe Roll Mode's Time Constant is %.2f\n",tau_r)
wn_r = 1./tau_r; %fprintf("\nThe Roll Mode's Natural Frequency is %.2f\n",wn_r)
wn_dr = sqrt(NB+(YB*Nr)/U1); %fprintf("\nThe Dutch Roll Mode's Natural Frequency is %.2f\n",wn_dr)
zeta_dr = -(Nr+YB/U1)/(2*wn_dr); %fprintf("\nThe Dutch Roll Mode's Damping Ratio is %.2f\n",zeta_dr)
tau_dr = 1/(wn_dr*zeta_dr); %fprintf("\nThe Dutch Roll Mode's Time Constant is %.2f\n",tau_dr)
zeta_s = 1; %fprintf("\nThe Spiral Mode's Damping Ratio is %.2f\n",zeta_s)
tau_s = -LB/(LB*Nr-NB*Lr); %fprintf("\nThe Spiral Mode's Time Constant is %.2f\n",tau_s)
wn_s = 1./tau_s; %fprintf("\nThe Sprial Mode's Natural Frequency is %.2f\n",wn_s)

% real values for modes:
poles_lat = eig(A); % finds poles (dutch roll is the first two, roll is the third row and spiral is last row)
zeta_lat_real = abs(cos(angle(poles_lat))); % matlab website is so goated ^.^
wn_lat_real = abs(real(poles_lat))./zeta_lat_real; % 1 and 2 = DR, 3 = roll, 4 = spiral
tau_lat_real = 1./(abs(real(poles_lat))); % same order as others

% problem 5
[b_lat,~] = ss2tf(A,B(:,1),C,D(:,1),1);
[b_lat(5:8,:),a_lat] = ss2tf(A,B,C,D,2); % 1 thru 4 is for aileron response, 5 thru 8 is for rudder response
% DB_DA_tf = tf(U1*b_lat(1,:),U1*a_lat)
DP_DA_tf = tf(U1*b_lat(2,:),U1*[a_lat,0]); % bank angle per aileron
% DR_DA_tf = tf(U1*b_lat(3,:),U1*a_lat)
% DPHI_DA_tf = tf(U1*b_lat(4,:),U1*a_lat)
% DB_DR_tf = tf(U1*b_lat(5,:),U1*a_lat)
% DP_DR_tf = tf(U1*b_lat(6,:),U1*a_lat)
DR_DR_tf = tf(U1*b_lat(7,:),U1*[a_lat,0]); % heading angle per rudder
% DPHI_DR_tf = tf(U1*b_lat(8,:),U1*a_lat)
% I commented out all the ones we didn't need so it wasnt cluttered

figure, rlocus(DP_DA_tf), title("Bank Angle to Aileron Transfer Function");
figure, rlocus(DR_DR_tf), title("Heading Angle to Rudder Transfer Function");
figure, bode(DP_DA_tf), margin(DP_DA_tf), title("Bank Angle to Aileron Transfer Function");
figure, bode(DR_DR_tf), margin(DR_DR_tf), title("Heading Angle to Rudder Transfer Function");

deflection = 2; % in degrees
t = 0:0.01:1000;
deltae = deg2rad(deflection)*ones(size(t));
[x,t]=lsim(sys_lon,deltae,t); 
figure, tiledlayout("Vertical"), nexttile
plot(t,x(:,1)+U1), ylabel("U, [ft/s]"), title("Response to Elevator Deflection \delta_e"), nexttile, 
plot(t,rad2deg(x(:,2)+AoA)), ylabel("\alpha, [deg]"), nexttile
plot(t,rad2deg(x(:,3))), ylabel("q, [deg/s]"), nexttile
plot(t,rad2deg(x(:,4))), ylabel("\theta, [deg]"), nexttile
plot(t, rad2deg(deltae)), ylabel("\delta_e, [deg]"), xlabel("Time(s)");

deltaa = deg2rad(30)*ones(size(t));
deltar = deg2rad(30)*ones(size(t));

% delta_a
[x,t]=lsim(sys_lat(:,1),deltaa,t); 
figure, tiledlayout("Vertical"), nexttile
plot(t,x(:,1)+U1), ylabel("\beta , [deg]"), title("Response to Aileron Deflection \delta_a"), nexttile, 
plot(t,rad2deg(x(:,2)+AoA)), ylabel("P, [deg/s]"), nexttile
plot(t,rad2deg(x(:,3))), ylabel("r, [deg/s]"), nexttile
plot(t,rad2deg(x(:,4)+AoA)), ylabel("\phi, [deg]"), nexttile
plot(t, rad2deg(deltaa)), ylabel("\delta_a, [deg]"), xlabel("Time(s)");

% delta_r
[x,t]=lsim(sys_lat(:,2),deltar,t); 
figure, tiledlayout("Vertical"), nexttile
plot(t,x(:,1)+U1), ylabel("\beta, [deg]"), title("Response to Rudder Deflection \delta_r"), nexttile, 
plot(t,rad2deg(x(:,2))), ylabel("P, [deg/s]"), nexttile
plot(t,rad2deg(x(:,3))), ylabel("r, [deg/s]"), nexttile
plot(t,rad2deg(x(:,4))), ylabel("\phi, [deg]"), nexttile
plot(t, rad2deg(deltar)), ylabel("\delta_r, [deg]"), xlabel("Time(s)");

fprintf("If you wanna see the Gain Margin and Phase Margin for a bode plot you need to delete the title for whatever plot, I prefer having titles so they don't get mixed up")