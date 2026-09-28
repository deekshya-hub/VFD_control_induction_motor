
% Dynamic Simulation and Constant V/f Control of a Three-Phase Induction Motor
% Reproduces the project results and generates four CV-quality figures.
%
% Model:
%   50 Hz, 338.8 V peak -> 40 Hz, 271.04 V peak at t = 5 s
%   Constant V/f ratio
%   20 N.m constant load torque
%
% The motor is represented in the synchronous dq reference frame.
% States:
%   psi_ds, psi_qs, psi_dr, psi_qr, omega_r

clear; clc; close all;

%% Motor parameters
p       = 2;          % pole pairs
Rs      = 1.77;       % stator resistance, ohm
Rr      = 1.34;       % rotor resistance, ohm
Lls     = 0.0139;     % stator leakage inductance, H
Llr     = 0.0121;     % rotor leakage inductance, H
Lm      = 0.3687;     % magnetizing inductance, H
J       = 0.001;      % inertia, kg.m^2
B       = 0;          % viscous damping, N.m.s/rad
Tload   = 20;         % load torque, N.m

Ls = Lls + Lm;
Lr = Llr + Lm;
Minv = inv([Ls Lm; Lm Lr]);

Vbase_pk = 338.8;     % phase-voltage peak at 50 Hz
fbase = 50;           % base frequency, Hz

%% Simulation settings
tspan = [0 10];
x0 = zeros(5,1);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-4);

[t,x] = ode45(@(t,x) motorODE(t,x,p,Rs,Rr,Ls,Lr,Lm,J,B,Tload,...
                   Vbase_pk,fbase,Minv), tspan, x0, opts);

%% Calculate outputs
speed_rpm = x(:,5)*60/(2*pi);
Te = zeros(size(t));
ids = zeros(size(t));
iqs = zeros(size(t));

for k = 1:length(t)
    psi_ds = x(k,1); psi_qs = x(k,2);
    psi_dr = x(k,3); psi_qr = x(k,4);

    i = Minv*[psi_ds; psi_dr];
    ids(k) = i(1);

    i = Minv*[psi_qs; psi_qr];
    iqs(k) = i(1);

    Te(k) = 1.5*p*(psi_ds*iqs(k)-psi_qs*ids(k));
end

%% Electrical angle and three phase currents
theta = zeros(size(t));
for k = 1:length(t)
    if t(k) < 5
        theta(k) = 2*pi*50*t(k);
    else
        theta(k) = 2*pi*50*5 + 2*pi*40*(t(k)-5);
    end
end

ia = ids.*cos(theta) - iqs.*sin(theta);
ib = ids.*cos(theta-2*pi/3) - iqs.*sin(theta-2*pi/3);
ic = ids.*cos(theta+2*pi/3) - iqs.*sin(theta+2*pi/3);

%% Results
fprintf('\n===== PROJECT RESULTS =====\n');
fprintf('Speed at 4.9 s: %.2f rpm\n', speed_rpm(find(t>=4.9,1)));
fprintf('Final speed at 40 Hz: %.2f rpm\n', speed_rpm(end));
fprintf('Final electromagnetic torque: %.2f N.m\n', mean(Te(t>9)));
fprintf('Final Phase-A RMS current: %.2f A\n', rms(ia(t>9)));
fprintf('V/f ratio: %.4f V/Hz (peak phase voltage)\n', Vbase_pk/fbase);

%% Plot 1: Speed
figure('Color','w');
plot(t,speed_rpm,'LineWidth',1.5);
grid on;
xline(5,'--','Frequency change');
xlabel('Time (s)');
ylabel('Motor speed (rpm)');
title('Induction Motor Speed Response under Constant V/f Control');
exportgraphics(gcf,'01_speed_response.png','Resolution',220);

%% Plot 2: Torque
figure('Color','w');
plot(t,Te,'LineWidth',1.5);
grid on;
xline(5,'--','Frequency change');
xlabel('Time (s)');
ylabel('Electromagnetic torque (N.m)');
title('Electromagnetic Torque Response');
exportgraphics(gcf,'02_torque_response.png','Resolution',220);

%% Plot 3: Three-phase current
idx = t >= 4.85 & t <= 5.15;
figure('Color','w');
plot(t(idx),ia(idx),'LineWidth',1.1); hold on;
plot(t(idx),ib(idx),'LineWidth',1.1);
plot(t(idx),ic(idx),'LineWidth',1.1);
xline(5,'--','Frequency change');
grid on;
xlabel('Time (s)');
ylabel('Phase current (A)');
title('Three-Phase Current During 50 Hz to 40 Hz Transition');
legend('I_a','I_b','I_c','Location','best');
exportgraphics(gcf,'03_three_phase_current.png','Resolution',220);

%% Plot 4: V/f characteristic
f = [50 40];
V = Vbase_pk*(f/fbase);
figure('Color','w');
plot(f,V,'o-','LineWidth',1.5);
grid on;
xlabel('Electrical frequency (Hz)');
ylabel('Phase-voltage peak (V)');
title('Constant V/f Characteristic');
exportgraphics(gcf,'04_vf_characteristic.png','Resolution',220);

%% Save numerical results
results = table(t,speed_rpm,Te,ia,ib,ic,...
    'VariableNames',{'Time_s','Speed_rpm','Torque_Nm','Ia_A','Ib_A','Ic_A'});
writetable(results,'simulation_results.csv');
save('simulation_results.mat','t','speed_rpm','Te','ia','ib','ic');

%% Local ODE function
function dx = motorODE(t,x,p,Rs,Rr,Ls,Lr,Lm,J,B,Tload,Vbase_pk,fbase,Minv)

psi_ds = x(1);
psi_qs = x(2);
psi_dr = x(3);
psi_qr = x(4);
omega_r = x(5);

if t < 5
    f = 50;
else
    f = 40;
end

omega_e = 2*pi*f;
Vpk = Vbase_pk*(f/fbase);

vds = 0;
vqs = Vpk;

i = Minv*[psi_ds;psi_dr];
ids = i(1);
idr = i(2);

i = Minv*[psi_qs;psi_qr];
iqs = i(1);
iqr = i(2);

slip_speed = omega_e - p*omega_r;

dpsi_ds = vds - Rs*ids + omega_e*psi_qs;
dpsi_qs = vqs - Rs*iqs - omega_e*psi_ds;
dpsi_dr = -Rr*idr + slip_speed*psi_qr;
dpsi_qr = -Rr*iqr - slip_speed*psi_dr;

Te = 1.5*p*(psi_ds*iqs - psi_qs*ids);
domega_r = (Te - Tload - B*omega_r)/J;

dx = [dpsi_ds;dpsi_qs;dpsi_dr;dpsi_qr;domega_r];
end
