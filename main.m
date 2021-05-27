clc
clear
close all
addpath('multicomplex')

%CW_trajectory([40 -400 500 0 0 1], 408e3, 5, 1000)
%cal_trajectory([-0.00018, -0.000002, 0.0012], [40 -400 500 0 0 1], 408e3, 5, 1000, 1e-15);
%[control_parameters,fval] = OptimiseTrajectory([-0.00018, -0.000002, 0.0012], [40 -400 500 0 0 1],  408e3, 5, 1000, 1e-15);
%BacktrackLineMethod([-0.00018, -0.000002, 0.0012], [40 -40 500 1 1 1],  408e3, 5, 1000, 1e-15);

%% Initializing Variables
%num_orbits = 0.82;
h = 1e-15;
dt = multicomplex(inputconverter(5,4,h)); % in seconds
num_steps = 1000; % round(P*num_orbits/dt)
animation_increments = 50;
reference_altitude = 408e3;
y0 = [40 -400 500 0 0 1]; % [x0 y0 z0 xdot0 ydot0 zdot0]

%RESULTING VARIABLES
mu = 3.986004418e14;
r0 = 6371e3 + reference_altitude;
omega = sqrt(mu/r0^3);
P = (2*pi*1/omega);
for k = 1:num_steps
    t_array(k) = k*dt;
end
%t_array = 0:dt:t_duration;

%% Calculate Exact Orbit
%{
ECEF = cylindrical_to_ECEF(y0,t(1),omega,r0); % Earth Centered Earth Fixed
keplerian_coord = ECEF_to_keplerian2D(ECEF,mu);

kepl_sol =  zeros(length(t_array),3);
for k = 1:size(t_array,2)
   kepl_sol(k,:) = keplerian_to_ECEF2D(keplerian_coord, mu, t_array(k)); % relative spacecraft to reference position
end
%}
%{
r = ECEF(1:3)./1000; % convert to km
v = ECEF(4:6)./1000; % convert to km/s
[a, ecc, incl, RAAN, argp, nu] = eci2orb_gooding(mu, r, v)
keplerian2ijk(a, ecc, incl, abs(RAAN), argp, nu)

h = cross(ECEF(1:3),ECEF(4:6));
nhat=cross([0 0 1],h);

evec = ((norm(v)^2-mu/norm(r))*r-dot(r,v)*v)/mu;
e = norm(evec);

energy = norm(v)^2/2-mu/norm(r);

if abs(e-1.0)>0.0001
   a = -mu/(2*energy);
   p = a*(1-e^2);
else
   p = norm(h)^2/mu;
   a = inf;
end

i = acos(h(3)/norm(h));

Omega = acos(n(1)/norm(n));
if n(2)<0
   Omega = 360-Omega;
end

argp = acos(dot(n,evec)/(norm(n)*e));

if e(3)<0
   argp = 360-argp;
end
nu = acos(dot(evec,r)/(e*norm(r)));

if dot(r,v)<0
   nu = 360 - nu;
end
[a, e, i, Omega, omega, nu]
%}

%% CW Solution (Clohessy-Wiltshire)
t_array_real = real(t_array);
CW_sol = zeros(length(t_array_real),6);
ECEF_spacecraft = zeros(length(t_array_real),3);
ECEF_reference = zeros(length(t_array_real),3);
for k = 1:size(t_array,2)
   CW_sol(k,:) = CW_solution(omega, t_array_real(k), y0); % relative spacecraft to reference position
   [ECEF_spacecraft(k,:),~] = cylindrical_to_ECEF(CW_sol(k,:),t_array_real(k),omega,r0);
   [ECEF_reference(k,:),~] = cylindrical_to_ECEF([0 0 0 0 0 0],t_array_real(k),omega,r0);
end

%% Control Solution
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*omega^2 0 0 0 2*omega  0;
     0 0 0 -2*omega 0 0;
     0 0 -omega^2 0 0 0];
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];
%{
x_spline = [-0.00019*multicomplex(inputconverter(1,1,1e-10)), -0.00018*multicomplex(inputconverter(1,2,1e-10))];
y_spline = [-0.000002*multicomplex(inputconverter(1,3,1e-10)), -0.000002*multicomplex(inputconverter(1,4,1e-10))];
z_spline = [0.0012*multicomplex(inputconverter(1,5,1e-10)), 0.0012*multicomplex(inputconverter(1,6,1e-10))];

syms t
x_control(t) = polyfit(linspace(0,t_array(end),length(x_spline)), x_spline, length(x_spline)-1);
y_control(t) = polyfit(linspace(0,t_array(end),length(y_spline)), y_spline, length(y_spline)-1);
z_control(t) = polyfit(linspace(0,t_array(end),length(z_spline)), z_spline, length(z_spline)-1);
%}

fd_state_traj.ax = zeros(length(t_array),6);
fd_state_traj.ay = zeros(length(t_array),6);
fd_state_traj.az = zeros(length(t_array),6);
fd_state_traj.dt = zeros(length(t_array),6);
fd_state_traj.axay = zeros(length(t_array),6);
fd_state_traj.axaz = zeros(length(t_array),6);
fd_state_traj.ayaz = zeros(length(t_array),6);
fd_state_traj.dtax = zeros(length(t_array),6);
fd_state_traj.dtay = zeros(length(t_array),6);
fd_state_traj.dtaz = zeros(length(t_array),6);
fd_state_traj.axn = zeros(length(t_array),6);
fd_state_traj.ayn = zeros(length(t_array),6);
fd_state_traj.azn = zeros(length(t_array),6);

fd_state_traj.ax(1,:) = y0;
fd_state_traj.ay(1,:) = y0;
fd_state_traj.az(1,:) = y0;
fd_state_traj.dt(1,:) = y0;
fd_state_traj.dtax(1,:) = y0;
fd_state_traj.dtay(1,:) = y0;
fd_state_traj.dtaz(1,:) = y0;
fd_state_traj.axay(1,:) = y0;
fd_state_traj.axaz(1,:) = y0;
fd_state_traj.ayaz(1,:) = y0;
fd_state_traj.axn(1,:) = y0;
fd_state_traj.ayn(1,:) = y0;
fd_state_traj.azn(1,:) = y0;

fd_control_traj.ax = zeros(length(t_array),3);
fd_control_traj.ay = zeros(length(t_array),3);
fd_control_traj.az = zeros(length(t_array),3);
fd_control_traj.axay = zeros(length(t_array),3);
fd_control_traj.axaz = zeros(length(t_array),3);
fd_control_traj.ayaz = zeros(length(t_array),3);
fd_control_traj.axn = zeros(length(t_array),3);
fd_control_traj.ayn = zeros(length(t_array),3);
fd_control_traj.azn = zeros(length(t_array),3);

fd_control_step = 1e-7;
fd_control_step_2 = 1e-7;


control_traj = cell(length(t_array),3);
control_traj(:,:) = {0};
complex_state_traj = cell(length(t_array),1);
complex_state_traj{1} = [multicomplex(y0(1));multicomplex(y0(2));multicomplex(y0(3));multicomplex(y0(4));multicomplex(y0(5));multicomplex(y0(6))];
real_state_traj = zeros(length(t_array),6);
real_state_traj(1,:) = y0;

control_perturbation.ax = multicomplex(inputconverter(0,1,h));
control_perturbation.ay = multicomplex(inputconverter(0,2,h));
control_perturbation.az = multicomplex(inputconverter(0,3,h));

for i = 1:round(length(t_array))
    control_traj{i,1} = -0.00018 + control_perturbation.ax;
    control_traj{i,2} = -0.000002 + control_perturbation.ay;
    control_traj{i,3} = 0.0012 + control_perturbation.az;
    
    fd_control_traj.ax(i,:) = [-0.00018+fd_control_step,-0.000002,0.0012];
    fd_control_traj.ay(i,:) = [-0.00018,-0.000002+fd_control_step,0.0012];
    fd_control_traj.az(i,:) = [-0.00018,-0.000002,0.0012+fd_control_step];
    fd_control_traj.dt(i,:) = [-0.00018,-0.000002,0.0012];
    fd_control_traj.axay(i,:) = [-0.00018+fd_control_step_2,-0.000002+fd_control_step_2,0.0012];
    fd_control_traj.axaz(i,:) = [-0.00018+fd_control_step_2,-0.000002,0.0012+fd_control_step_2];
    fd_control_traj.ayaz(i,:) = [-0.00018,-0.000002+fd_control_step_2,0.0012+fd_control_step_2];
    fd_control_traj.axn(i,:) = [-0.00018-fd_control_step,-0.000002,0.0012];
    fd_control_traj.ayn(i,:) = [-0.00018,-0.000002-fd_control_step,0.0012];
    fd_control_traj.azn(i,:) = [-0.00018,-0.000002,0.0012-fd_control_step];
end

STM = (eye(6) + dt*A + dt^2*A^2/2 + dt^3*A^3/6 + dt^4*A^4/24); % State Transition Matrix
CTM = (dt*eye(6) + dt^2*A/2 + dt^3*A^2/6 + dt^4*A^3/24); % Control Transition Matrix inv(A)*(STM - eye(6))
STM_real = real(STM);
CTM_real = real(CTM);
dt_fd = real(dt) + fd_control_step;
STM_fd = (eye(6) + dt_fd*A + dt_fd^2*A^2/2 + dt_fd^3*A^3/6 + dt_fd^4*A^4/24); % State Transition Matrix
CTM_fd = (dt_fd*eye(6) + dt_fd^2*A/2 + dt_fd^3*A^2/6 + dt_fd^4*A^3/24); % Control Transition Matrix inv(A)*(STM - eye(6))

for k = 2:length(t_array)
    complex_state_traj{k} = STM*complex_state_traj{k-1} + CTM*B*[control_traj{k,:}]';
    real_state_traj(k,:) = real(complex_state_traj{k});
    
    fd_state_traj.ax(k,:) = STM_real*fd_state_traj.ax(k-1,:)' + CTM_real*B*fd_control_traj.ax(k,:)';
    fd_state_traj.ay(k,:) = STM_real*fd_state_traj.ay(k-1,:)' + CTM_real*B*fd_control_traj.ay(k,:)';
    fd_state_traj.az(k,:) = STM_real*fd_state_traj.az(k-1,:)' + CTM_real*B*fd_control_traj.az(k,:)';
    fd_state_traj.dt(k,:) = STM_fd*fd_state_traj.dt(k-1,:)' + CTM_fd*B*fd_control_traj.dt(k,:)'; %d/d_dt
    fd_state_traj.dtax(k,:) = STM_fd*fd_state_traj.dtax(k-1,:)' + CTM_fd*B*fd_control_traj.ax(k,:)'; %d^2/d_(dt_ax)
    fd_state_traj.dtay(k,:) = STM_fd*fd_state_traj.dtay(k-1,:)' + CTM_fd*B*fd_control_traj.ay(k,:)'; %d^2/d_(dt_ay)
    fd_state_traj.dtaz(k,:) = STM_fd*fd_state_traj.dtaz(k-1,:)' + CTM_fd*B*fd_control_traj.az(k,:)'; %d^2/d_(dt_az)
    fd_state_traj.axay(k,:) = STM_fd*fd_state_traj.axay(k-1,:)' + CTM_fd*B*fd_control_traj.axay(k,:)'; %d^2/d_(ax_ay)
    fd_state_traj.axaz(k,:) = STM_fd*fd_state_traj.axaz(k-1,:)' + CTM_fd*B*fd_control_traj.axaz(k,:)'; %d^2/d_(ax%_az)
    fd_state_traj.ayaz(k,:) = STM_fd*fd_state_traj.ayaz(k-1,:)' + CTM_fd*B*fd_control_traj.ayaz(k,:)'; %d^2/d_(ay_az)
    fd_state_traj.axn(k,:) = STM_real*fd_state_traj.axn(k-1,:)' + CTM_real*B*fd_control_traj.axn(k,:)';
    fd_state_traj.ayn(k,:) = STM_real*fd_state_traj.ayn(k-1,:)' + CTM_real*B*fd_control_traj.ayn(k,:)';
    fd_state_traj.azn(k,:) = STM_real*fd_state_traj.azn(k-1,:)' + CTM_real*B*fd_control_traj.azn(k,:)';
end

%% Proportional Controler
%{
control_gain = 0.000001;
proportional_traj = zeros(length(t_array),6);
proportional_traj(1,:) = y0;
proportional_controler = zeros(length(t_array),3);

for k = 2:length(t_array)
    proportional_controler(k,:) = [-proportional_traj(k-1,1),proportional_traj(k-1,2),proportional_traj(k-1,2)].*control_gain;
    proportional_traj(k,:) = STM*proportional_traj(k-1,:)' + CTM*B*proportional_controler(k,:)';
end
proportional_cost = sum_control_cost(proportional_traj,dt);
%}

%% Euler and Runge-Kutta Methods
%{
euler_sol = zeros(length(t_array),6);
euler_sol(1,:) = y0;
RK_sol = zeros(length(t_array),6);
RK_sol(1,:) = y0;
tay_sol = zeros(length(t_array),6);
tay_sol(1,:) = y0;

for k = 2:length(t_array)
    %Euler
    euler_sol(k,:) = ( (eye(6) + dt*A)*euler_sol(k-1,:)' )';
    %Runge-Kutta
    k1 = dt*A*RK_sol(k-1,:)';
    k2 = dt*A*(RK_sol(k-1,:)' + 0.5*k1);
    k3 = dt*A*(RK_sol(k-1,:)' + 0.5*k2);
    k4 = dt*A*(RK_sol(k-1,:)' + k3);
    RK_sol(k,:) = ( RK_sol(k-1,:)' + 1/6*(k1 + 2*k2 + 2*k3 + k4) )';
    % Fourth Order Taylor expansion
    tay_sol(k,:) = ( (eye(6) + dt*A + dt^2*A^2/2 + dt^3*A^3/6 + dt^4*A^4/24)*tay_sol(k-1,:)' )';
end
%}

%% Final State and Cost Functions
final_state = strings(size(complex_state_traj{end},1),1);
for i = 1:size(complex_state_traj{end},1)
    final_state(i) = repr(complex_state_traj{end}(i));
end

complex_sensitivity = zeros(6,7,length(complex_state_traj));

for k = 1:length(complex_state_traj)
    for i = 1:6
        complex_sensitivity(i,1,k) = CXn(complex_state_traj{k}(i),1)/CXn(control_perturbation.ax,1);
        complex_sensitivity(i,2,k) = CXn(complex_state_traj{k}(i),2)/CXn(control_perturbation.ay,2);
        complex_sensitivity(i,3,k) = CXn(complex_state_traj{k}(i),3)/CXn(control_perturbation.az,3);
        complex_sensitivity(i,4,k) = CXn(complex_state_traj{k}(i),4)/CXn(dt,4);

        complex_sensitivity(i,5,k) = CXn(complex_state_traj{k}(i),[1 4])/CXn(dt,4)^2; % 'd^2/d_dt_ax'
        complex_sensitivity(i,6,k) = CXn(complex_state_traj{k}(i),[2 4])/CXn(dt,4)^2; % 'd^2/d_dt_ay'
        complex_sensitivity(i,7,k) = CXn(complex_state_traj{k}(i),[3 4])/CXn(dt,4)^2; % 'd^2/d_dt_az'
    end
end

final_complex_sensitivity = array2table(complex_sensitivity(:,:,end),'VariableNames',{'d/d_ax','d/d_ay','d/d_az','d/d_dt','d^2/d_(dt_ax)','d^2/d_(dt_ay)','d^2/d_(dt_az)'});
final_complex_sensitivity = [table({'x';'y';'z';'x_dot';'y_dot';'z_dot'},'VariableNames',{' '}), final_complex_sensitivity];


fd_sensitivity(:,1) = (fd_state_traj.ax(end,:) - real_state_traj(end,:))'./fd_control_step;
fd_sensitivity(:,2) = (fd_state_traj.ay(end,:) - real_state_traj(end,:))'./fd_control_step;
fd_sensitivity(:,3) = (fd_state_traj.az(end,:) - real_state_traj(end,:))'./fd_control_step;
fd_sensitivity(:,4) = (fd_state_traj.dt(end,:) - real_state_traj(end,:))'./fd_control_step;
fd_sensitivity(:,5) = (fd_state_traj.dtax(end,:) - fd_state_traj.dt(end,:) - fd_state_traj.ax(end,:) + real_state_traj(end,:))'./fd_control_step^2;
fd_sensitivity(:,6) = (fd_state_traj.dtay(end,:) - fd_state_traj.dt(end,:) - fd_state_traj.ay(end,:) + real_state_traj(end,:))'./fd_control_step^2;
fd_sensitivity(:,7) = (fd_state_traj.dtaz(end,:) - fd_state_traj.dt(end,:) - fd_state_traj.az(end,:) + real_state_traj(end,:))'./fd_control_step^2;

fd_sensitivity = array2table(fd_sensitivity,'VariableNames',{'d/d_ax','d/d_ay','d/d_az','d/d_dt','d^2/d_(dt_ax)','d^2/d_(dt_ay)','d^2/d_(dt_az)'});
fd_sensitivity = [table({'x';'y';'z';'x_dot';'y_dot';'z_dot'},'VariableNames',{' '}), fd_sensitivity];


control_cost = sum_control_cost(control_traj,dt);
traj_duration = dt*num_steps;


%% Objective Gradients
% Complex Method Gradient
complex_radius = multicomplex(complex_state_traj{end}(1).zn.^2 + complex_state_traj{end}(2).zn.^2 + complex_state_traj{end}(3).zn.^2);
objective_value = real(complex_radius);
complex_obj_sensitivity = [CXn(complex_radius,1)/CXn(control_perturbation.ax,1);CXn(complex_radius,2)/CXn(control_perturbation.ay,2);CXn(complex_radius,3)/CXn(control_perturbation.az,3)];

% Complex Method Hessian
complex_hessian = zeros(3);
complex_hessian(1,2) = CXn(complex_radius,[1 2])/CXn(control_perturbation.ax,1)^2;
complex_hessian(1,3) = CXn(complex_radius,[1 3])/CXn(control_perturbation.ay,2)^2;
complex_hessian(2,3) = CXn(complex_radius,[2 3])/CXn(control_perturbation.az,3)^2;

% Finite Difference Gradient
real_radius = real_state_traj(end,1)^2 + real_state_traj(end,2)^2 + real_state_traj(end,3)^2;
fd_radius.ax = fd_state_traj.ax(end,1)^2 + fd_state_traj.ax(end,2)^2 + fd_state_traj.ax(end,3)^2;
fd_radius.ay = fd_state_traj.ay(end,1)^2 + fd_state_traj.ay(end,2)^2 + fd_state_traj.ay(end,3)^2;
fd_radius.az = fd_state_traj.az(end,1)^2 + fd_state_traj.az(end,2)^2 + fd_state_traj.az(end,3)^2;

fd_obj_sensitivity = [(fd_radius.ax - real_radius)/fd_control_step ; (fd_radius.ay - real_radius)/fd_control_step; (fd_radius.az - real_radius)/fd_control_step];

% Finite Difference Hessian
fd_radius.axay = fd_state_traj.axay(end,1)^2 + fd_state_traj.axay(end,2)^2 + fd_state_traj.axay(end,3)^2;
fd_radius.axaz = fd_state_traj.axaz(end,1)^2 + fd_state_traj.axaz(end,2)^2 + fd_state_traj.axaz(end,3)^2;
fd_radius.ayaz = fd_state_traj.ayaz(end,1)^2 + fd_state_traj.ayaz(end,2)^2 + fd_state_traj.ayaz(end,3)^2;
fd_radius.axn = fd_state_traj.axn(end,1)^2 + fd_state_traj.axn(end,2)^2 + fd_state_traj.axn(end,3)^2;
fd_radius.ayn = fd_state_traj.ayn(end,1)^2 + fd_state_traj.ayn(end,2)^2 + fd_state_traj.ayn(end,3)^2;
fd_radius.azn = fd_state_traj.azn(end,1)^2 + fd_state_traj.azn(end,2)^2 + fd_state_traj.azn(end,3)^2;

fd_obj_hessian = zeros(3);

fd_obj_hessian(1,1) = (fd_radius.ax - 2*real_radius + fd_radius.axn)/fd_control_step_2^2;
fd_obj_hessian(2,2) = (fd_radius.ay - 2*real_radius + fd_radius.ayn)/fd_control_step_2^2;
fd_obj_hessian(3,3) = (fd_radius.az - 2*real_radius + fd_radius.azn)/fd_control_step_2^2;

fd_obj_hessian(1,2) = (fd_radius.axay - fd_radius.ax - fd_radius.ay + real_radius)/fd_control_step_2^2;
fd_obj_hessian(2,1) = fd_obj_hessian(1,2);
fd_obj_hessian(1,3) = (fd_radius.axaz - fd_radius.ax - fd_radius.az + real_radius)/fd_control_step_2^2;
fd_obj_hessian(3,1) = fd_obj_hessian(1,3);
fd_obj_hessian(2,3) = (fd_radius.ayaz - fd_radius.ay - fd_radius.az + real_radius)/fd_control_step_2^2;
fd_obj_hessian(3,2) = fd_obj_hessian(2,3);



%% Plotting Inanimate Curves
figure('Renderer', 'painters', 'Position', [400 10 1300 900])
hold on
t = tiledlayout(2,2);
t.TileSpacing = 'none';
t.Padding = 'tight';

nexttile(2)
hold on
plot3(CW_sol(:,2),CW_sol(:,1),CW_sol(:,3),'r');
plot3(real_state_traj(:,2),real_state_traj(:,1),real_state_traj(:,3),'b');
plot3(0,0,0,'og','LineWidth',2)
%plot(proportional_traj(:,2),proportional_traj(:,1),'k')

nexttile(4)
hold on
plot3(CW_sol(:,5),CW_sol(:,4),CW_sol(:,6),'r');
plot3(real_state_traj(:,5),real_state_traj(:,4),real_state_traj(:,6),'b');
plot3(0,0,0,'og','LineWidth',2)
%plot(proportional_traj(:,5),proportional_traj(:,4),'k')

%{
% PLOT ERROR BETWEEN EULER AND CW AND RK
figure(3)
hold on
plot(100*(CW_sol(:,1)-euler_sol(:,1))./CW_sol(:,1));
plot(100*(CW_sol(:,2)-euler_sol(:,2))./CW_sol(:,2));
plot(100*(RK_sol(:,1)-RK_sol(:,1))./RK_sol(:,1));
plot(100*(RK_sol(:,2)-RK_sol(:,2))./RK_sol(:,2));
legend('euler y error','euler x error','RK y error','RK x error');
xlabel('time (s)')
ylabel('% error')
figure(1)
%}
%{
%% PLOT ERROR BETWEEN KEPLER SOLUTION AND PROPAGATED SOLUTION
figure(4)
hold on
plot(earth_spacecraft_orbit(:,1));
plot(kepl_sol(:,1));
plot(earth_spacecraft_orbit(:,2));
plot(kepl_sol(:,2));
x_error = 100*(earth_spacecraft_orbit(:,1)-kepl_sol(:,1))./earth_spacecraft_orbit(:,1);
y_error = 100*(earth_spacecraft_orbit(:,2)-kepl_sol(:,2))./earth_spacecraft_orbit(:,2);
plot(x_error)
plot(y_error)
figure(3)
legend('CW x coordinate', 'kepler x coordinate')
xlabel('time (s)')
ylabel('position m')
%}

figure(2)
hold on
title('XY Sensitivity Quiver Plot')
plot(real_state_traj(:,2),real_state_traj(:,1),'b');
quiver_vectors = [real_state_traj(:,2),real_state_traj(:,1),reshape(complex_sensitivity(2,1,:),[],1),reshape(complex_sensitivity(1,1,:),[],1)];
quiver_vectors = quiver_vectors(1:10:size(quiver_vectors,1),:);
quiver(quiver_vectors(:,1),quiver_vectors(:,2),quiver_vectors(:,3),quiver_vectors(:,4)*10,0.1);
ylabel('$x_{rel}$ (m)','Interpreter','latex')
xlabel('$y_{rel}$ (m)','Interpreter','latex')
set(gca, 'xdir', 'reverse')
figure(1)


%% Plot z state space transition plot
Z_matrix = [0       1;
           -omega^2 0];

z_dot_scale = 1000;
border = 0.5;
z_range = range([real_state_traj(:,3);CW_sol(:,3)]);
z_dot_range = range([real_state_traj(:,6);CW_sol(:,6)]);

min_z = min([real_state_traj(:,3);CW_sol(:,3)])-border*z_range;
min_z_dot = min([real_state_traj(:,6);CW_sol(:,6)])-border*z_dot_range;
max_z = max([real_state_traj(:,3);CW_sol(:,3)])+border*z_range;
max_z_dot = max([real_state_traj(:,6);CW_sol(:,6)])+border*z_dot_range;

z_intervals = (max_z-min_z)/18;
z_dot_intervals = (max_z_dot-min_z_dot)/18;

x_range = min_z:z_intervals:max_z;
y_range = max_z_dot:-z_dot_intervals:min_z_dot; % y_range = 20:-2:-20;
U = zeros(length(y_range),length(x_range));
V = zeros(length(y_range),length(x_range));
[z,z_dot] = meshgrid(x_range,y_range);

for i = 1:size(z,1)
    for j = 1:size(z,2)
        U(i,j) = Z_matrix(1,:)*[z(i,j);z_dot(i,j)];
        V(i,j) = Z_matrix(2,:)*[z(i,j);z_dot(i,j)];
    end
end

nexttile(1)
hold on
quiver(z,z_dot*z_dot_scale,U,V*z_dot_scale,'color',[0.25, 0.25, 0.25]) % quiver(z,z_dot*100,U,V*10000) 
plot(CW_sol(:,3),CW_sol(:,6)*z_dot_scale,'r') % plot(CW_sol(:,3),CW_sol(:,6)*1000,'r')
plot(real_state_traj(:,3),real_state_traj(:,6)*1000,'b') % plot(real_state_traj(:,3),real_state_traj(:,6)*1000,'b')
plot3(0,0,0,'og','LineWidth',2)
xlim([min_z max_z])
ylim([min_z_dot*z_dot_scale max_z_dot*z_dot_scale])
yt = get(gca,'ytick');    
set(gca,'YTick',yt, 'yticklabel',yt/1000)

%plot(proportional_traj(:,3),proportional_traj(:,6)*z_dot_scale,'k')

%% Plot x y state space transition plots
%{
XY_matrix = [0         0  1       0;
             0         0  0       1;
             3*omega^2 0  0       2*omega;
             0         0 -2*omega 0];

%z_dot_scale = 1000;
%border = 0.5;
y_range = range([real_state_traj(:,2);CW_sol(:,2)]);
x_range = range([real_state_traj(:,1);CW_sol(:,1)]);

min_y = min([real_state_traj(:,2);CW_sol(:,2)])-border*y_range;
min_x = min([real_state_traj(:,1);CW_sol(:,1)])-border*x_range;
max_y = max([real_state_traj(:,2);CW_sol(:,2)])+border*y_range;
max_x = max([real_state_traj(:,1);CW_sol(:,1)])+border*x_range;

y_intervals = (max_y-min_y)/18;
x_intervals = (max_x-min_x)/18;

x_range = min_y:y_intervals:max_y;
y_range = max_x:-x_intervals:min_x; % y_range = 20:-2:-20;
U = zeros(length(x_range),length(y_range));
V = zeros(length(x_range),length(y_range));
[y,x] = meshgrid(y_range,x_range);

for i = 1:size(z,1)
    for j = 1:size(z,2)
        U(i,j) = XY_matrix(1,:)*[x(i,j);y(i,j);x(i,j);y(i,j)];
        V(i,j) = XY_matrix(2,:)*[x(i,j);y(i,j);x(i,j);y(i,j)];
    end
end
%{
nexttile(1)
hold on
quiver(z,z_dot*z_dot_scale,U,V*z_dot_scale,'color',[0.25, 0.25, 0.25]) % quiver(z,z_dot*100,U,V*10000) 
plot(CW_sol(:,3),CW_sol(:,6)*z_dot_scale,'r') % plot(CW_sol(:,3),CW_sol(:,6)*1000,'r')
plot(real_state_traj(:,3),real_state_traj(:,6)*1000,'b') % plot(real_state_traj(:,3),real_state_traj(:,6)*1000,'b')
plot3(0,0,0,'og','LineWidth',2)
xlim([min_z max_z])
ylim([min_z_dot*z_dot_scale max_z_dot*z_dot_scale])
yt = get(gca,'ytick');    
set(gca,'YTick',yt, 'yticklabel',yt/1000)
%}
%}

%% Ploting sphere and orbit
nexttile(3)
hold on
sphere_size = 6371e3;
[X,Y,Z] = sphere(14);
surf(X*sphere_size,Y*sphere_size,Z*sphere_size);
%plot3(ECEF_reference(:,1),ECEF_reference(:,2),ECEF_reference(:,3),'g');
plot3(ECEF_spacecraft(:,1),ECEF_spacecraft(:,2),ECEF_spacecraft(:,3),'b');

title('Orbiting Reference')
axis equal
view([20 75])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('$z$ (m)','Interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0)

%% Plot Settings
nexttile(1)
xlabel('$z_{rel}$ (m)','Interpreter','latex')
ylabel('$\dot z_{rel}$ (m)','Interpreter','latex')
%set(get(gca,'YLabel'),'Rotation',0)
title('Relative Z Axis State Space')
legend('Zero Control Vector Field','Zero Control Trajectory', 'Controlled Trajectory','Target','AutoUpdate','off') 

nexttile(2)
set(gca, 'xdir', 'reverse')
title('Relative Position Trajectory')
ylabel('$x_{rel}$ (m)','Interpreter','latex')
xlabel('$y_{rel}$ (m)','Interpreter','latex')
zlabel('$z_{rel}$ (m)','Interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0)
%x_border = range(real_state_traj(:,2))*0.1 + 1;
%y_border = range(real_state_traj(:,1))*0.1 + 1;
%z_border = range(real_state_traj(:,3))*0.1 + 1;
%xlim([min(state_traj(:,2))-x_border max(state_traj(:,2))+x_border])
%ylim([min(state_traj(:,1))-y_border max(state_traj(:,1))+y_border])
%zlim([min(state_traj(:,3))-z_border max(state_traj(:,3))+z_border])
%view([-10 65])
legend('Zero Control Trajectory', 'Controlled Trajectory','Target','AutoUpdate','off') 

nexttile(4)
title('Relative Velocity Trajectory')
%view([-10 65])
%axis equal
ylabel('$\dot x_{rel}$ (m)','Interpreter','latex')
xlabel('$\dot y_{rel}$ (m)','Interpreter','latex')
zlabel('$\dot z_{rel}$ (m)','Interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0)
set(gca, 'xdir', 'reverse')
x_border = range(real_state_traj(:,4))*0.1 + 0.5;
y_border = range(real_state_traj(:,5))*0.1 + 0.5;
%xlim([min(state_traj(:,4))-x_border max(state_traj(:,4))+x_border])
%ylim([min(state_traj(:,5))-y_border max(CW_sol(:,5))+y_border])
legend('Zero Control Trajectory', 'Controlled Trajectory','Target','AutoUpdate','off') 

%% Animations
figure(1)
nexttile(1)
animated_z_no_control = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_z_with_control = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);

nexttile(2)
animated_relative_pos = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_relative_pos_traj = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);

nexttile(3)
animated_earth_reference_dot = animatedline('Marker','o','MarkerSize',9,'Color','g', 'MaximumNumPoints',1);
animated_earth_spacecraft_dot = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);

nexttile(4)
animated_relative_speed = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_relative_speed_traj = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);

%{
nexttile(5)
animated_earth_spacecraft_orbit = animatedline('LineWidth',1,'Color','b', 'MaximumNumPoints',2);
animated_earth_reference_orbit = animatedline('LineWidth',1,'Color','g', 'MaximumNumPoints',2);
animated_earth_spacecraft_dot2 = animatedline('Marker','x','MarkerSize',4,'LineWidth',2,'Color','b', 'MaximumNumPoints',1);
animated_earth_reference_dot2 = animatedline('Marker','o','MarkerSize',4,'LineWidth',2,'Color','g', 'MaximumNumPoints',1);
%}

for k = 2:animation_increments:length(t_array)
    nexttile(1)
    addpoints(animated_z_no_control,CW_sol(k,3),CW_sol(k,6)*z_dot_scale);
    addpoints(animated_z_with_control,real_state_traj(k,3),real_state_traj(k,6)*z_dot_scale);

    nexttile(2)
    addpoints(animated_relative_pos,CW_sol(k,2),CW_sol(k,1),CW_sol(k,3));
    addpoints(animated_relative_pos_traj,real_state_traj(k,2),real_state_traj(k,1),real_state_traj(k,3));
   
    nexttile(3)
    addpoints(animated_earth_reference_dot,ECEF_reference(k,1),ECEF_reference(k,2),ECEF_reference(k,3));
    addpoints(animated_earth_spacecraft_dot,ECEF_spacecraft(k,1),ECEF_spacecraft(k,2),ECEF_spacecraft(k,3));
    
    nexttile(4)
    addpoints(animated_relative_speed,CW_sol(k,5),CW_sol(k,4),CW_sol(k,6));
    addpoints(animated_relative_speed_traj,real_state_traj(k,5),real_state_traj(k,4),real_state_traj(k,6));
    
    %{
    nexttile(5)
    addpoints(animated_earth_reference_orbit,ECEF_reference(k,1),ECEF_reference(k,2),ECEF_reference(k,3));
    addpoints(animated_earth_spacecraft_orbit,ECEF_spacecraft(k,1),ECEF_spacecraft(k,2),ECEF_spacecraft(k,3));

    
    refX_mid = ECEF_reference(k,1) - 0.01*(ECEF_reference(k,1) - ECEF_reference(k-1,1));
    refY_mid = ECEF_reference(k,2) - 0.01*(ECEF_reference(k,2) - ECEF_reference(k-1,2));
    craftX_mid = ECEF_spacecraft(k,1) - 0.01*(ECEF_spacecraft(k,1) - ECEF_spacecraft(k-1,1));
    craftY_mid = ECEF_spacecraft(k,2) - 0.01*(ECEF_spacecraft(k,2) - ECEF_spacecraft(k-1,2));
    
    addpoints(animated_earth_reference_orbit,refX_mid,refY_mid,ECEF_reference(k,3));
    addpoints(animated_earth_spacecraft_orbit,craftX_mid,craftY_mid,ECEF_spacecraft(k,3));
    
    addpoints(animated_earth_reference_dot2,ECEF_reference(k,1),ECEF_reference(k,2),ECEF_reference(k,3));
    addpoints(animated_earth_spacecraft_dot2,ECEF_spacecraft(k,1),ECEF_spacecraft(k,2),ECEF_spacecraft(k,3));
    %}
    drawnow
end


t_array = t_array';


%% Functions
%{
function dydt = odefun(t,y,omega)
    dydt = zeros(6,1);
    dydt(1) = y(4);
    dydt(2) = y(5);
    dydt(3) = y(6);
    dydt(4) = 3*omega^2*y(1) + 2*omega*y(5);
    dydt(5) = -2*omega*y(5);
    dydt(6) = -omega^2*y(3);
end
%}

%initial_control_parameters = [-0.00018, -0.000002, 0.0012];

function [control_parameters,fval] = OptimiseTrajectory(initial_control_parameters,t0_state, ref_altitude, intial_time_step, num_steps, h)
    options = optimoptions('fmincon','SpecifyObjectiveGradient',false,'Display','iter'); % 'OutputFcn'?
    [control_parameters,fval] = fmincon(@(control_parameters)cal_trajectory(control_parameters, t0_state, ref_altitude, intial_time_step, num_steps, h),initial_control_parameters,[],[],[],[],[-.001;-.001;-.01],[.001;.001;.01],[],options);
end

function BacktrackLineMethod(initial_control_parameters,t0_state, ref_altitude, intial_time_step, num_steps, h)
    uiFig = uifigure('Position',[100 100 1125 475]);
    uT = uitable(uiFig);
    uT.ColumnName ={'iter', 'fun calls', 'step size t', 'objective val','x_f', 'y_f', 'z_f', 'ax', 'ay', 'az',};
    uT.Units='Normalized';
    uT.Position = [.1 .1 .8 .8];
    drawnow
    
    control_parameters = initial_control_parameters;
    beta = 0.3;
    alpha = 0.5;
    [objective, gradient,final_state] = calc_obj_grad(control_parameters);
    uT.Data(1,:) = [0 1 0 objective final_state(1:3) control_parameters(1:3)];
    drawnow
    for i = 1:10
        t = 1e-14;
        f_calls = 1;
                
        % backtracking line search
        while true
        [offset_objective, offset_gradient, final_state] = calc_obj_grad(control_parameters-t*gradient);
            if offset_objective > objective - sum(gradient.^2)*t*alpha
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        objective = offset_objective;
        control_parameters = control_parameters - t*gradient; %2e-16
        gradient = offset_gradient;
        uT.Data(1+i,:) = [i f_calls t objective final_state(1:3) control_parameters(1:3)];
        drawnow
    end
    
    function [obj, grad, state] = calc_obj_grad(control_parameters)
        [obj, grad, state] = cal_trajectory(control_parameters, t0_state, ref_altitude, intial_time_step, num_steps, h);
    end

end

function [objective_value, fun_gradient, final_state] = cal_trajectory(control_parameters, t0_state, ref_altitude, time_step, num_steps, h)
    dt = multicomplex(inputconverter(time_step,4,h));
    
    mu = 3.986004418e14;
    r0 = 6371e3 + ref_altitude;
    omega = sqrt(mu/r0^3);
    P = (2*pi*1/omega);
    for k = 1:num_steps
        t_array(k) = k*dt;
    end
    
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         3*omega^2 0 0 0 2*omega  0;
         0 0 0 -2*omega 0 0;
         0 0 -omega^2 0 0 0];
    B = [0 0 0;
         0 0 0;
         0 0 0;
         1 0 0;
         0 1 0;
         0 0 1];
     
    control_traj = cell(length(t_array),3);
    control_traj(:,:) = {0};
    complex_state_traj = cell(length(t_array),1);
    complex_state_traj{1} = [multicomplex(t0_state(1));multicomplex(t0_state(2));multicomplex(t0_state(3));multicomplex(t0_state(4));multicomplex(t0_state(5));multicomplex(t0_state(6))];
    real_state_traj = zeros(length(t_array),6);
    real_state_traj(1,:) = t0_state;
    control_perturbation.ax = multicomplex(inputconverter(0,1,h));
    control_perturbation.ay = multicomplex(inputconverter(0,2,h));
    control_perturbation.az = multicomplex(inputconverter(0,3,h));
    
    for i = 1:round(length(t_array))
        control_traj{i,1} = control_parameters(1) + control_perturbation.ax;
        control_traj{i,2} = control_parameters(2) + control_perturbation.ay;
        control_traj{i,3} = control_parameters(3) + control_perturbation.az;
    end
    
    STM = (eye(6) + dt*A + dt^2*A^2/2 + dt^3*A^3/6 + dt^4*A^4/24); % State Transition Matrix
    CTM = (dt*eye(6) + dt^2*A/2 + dt^3*A^2/6 + dt^4*A^3/24); % Control Transition Matrix inv(A)*(STM - eye(6))
    STM_real = real(STM);
    CTM_real = real(CTM);
    
    
    for k = 2:length(t_array)
        complex_state_traj{k} = STM*complex_state_traj{k-1} + CTM*B*[control_traj{k,:}]';
        real_state_traj(k,:) = real(complex_state_traj{k});
    end
    
    %final_pos_sum = multicomplex(abs(complex_state_traj{end}(1).zn) + abs(complex_state_traj{end}(2).zn) + abs(complex_state_traj{end}(3).zn));
    %objective_value = real(final_pos_sum);
    %final_pos_sum = complex_state_traj{end}(1) + complex_state_traj{end}(2) + complex_state_traj{end}(3);
    %objective_value = abs(real(complex_state_traj{end}(1))) + abs(real(complex_state_traj{end}(2))) + abs(real(complex_state_traj{end}(3)));
    complex_radius = complex_state_traj{end}(1)^2 + complex_state_traj{end}(2)^2 + complex_state_traj{end}(3)^2;
    objective_value = real(complex_radius);
    fun_gradient = [CXn(complex_radius,1)/CXn(control_perturbation.ax,1),CXn(complex_radius,2)/CXn(control_perturbation.ay,2),CXn(complex_radius,3)/CXn(control_perturbation.az,3)];

    %disp(control_parameters)
    final_state = real_state_traj(end,:);
    %disp(final_state)
end

function output = CW_trajectory(t0_state, ref_altitude, time_step, num_steps)
    dt = time_step;
    
    mu = 3.986004418e14;
    r0 = 6371e3 + ref_altitude;
    omega = sqrt(mu/r0^3);
    P = (2*pi*1/omega);
    for k = 1:num_steps
        t_array(k) = k*dt;
    end
    
    CW_sol = zeros(length(t_array),6);
    for k=1:length(t_array)
        CW_sol(k,:) = CW_solution(omega, t_array(k), t0_state);
    end
end

function output = CW_solution(omega,t, X0)
    x0 = X0(1);
    y0 = X0(2);
    z0 = X0(3);
    xdot0 = X0(4);
    ydot0 = X0(5);
    zdot0 = X0(6);
      
       
    xy_mat = [4-3*cos(omega*t), 0, sin(omega*t)/omega, 2*(1-cos(omega*t))/omega;
        6*(sin(omega*t)-omega*t),1, 2*(cos(omega*t)-1)/omega, (4*sin(omega*t)-3*omega*t)/omega;
        3*omega*sin(omega*t), 0, cos(omega*t), 2*sin(omega*t);
        6*omega*(cos(omega*t)-1), 0, -2*sin(omega*t), 4*cos(omega*t)-3];
    z_mat = [cos(omega*t), sin(omega*t)/omega;
        -omega*sin(omega*t), cos(omega*t)];
    
    local_sol = [xy_mat*[x0;y0;xdot0;ydot0]; z_mat*[z0;zdot0]];
    x = local_sol(1);
    y = local_sol(2);
    z = local_sol(5);
    xdot = local_sol(3);
    ydot = local_sol(4);
    zdot = local_sol(6);
    output = [x,y,z,xdot,ydot,zdot];
end

function [pos,velocity] = cylindrical_to_ECEF(rel_coord,t,omega,r0)
    nu = t*omega;
    rel_cylindrical = [rel_coord, r0]; % [x y z x_dot y_dot z_dot r0]
      
    ECEF = [cos(nu),       -sin(nu),       0, 0,       0,        0, cos(nu);
           sin(nu),        cos(nu),        0, 0,       0,        0, sin(nu);
           0,              0,              1, 0,       0,        0, 0;
           -omega*sin(nu), -omega*cos(nu), 0, cos(nu), -sin(nu), 0, -omega*sin(nu);
           omega*cos(nu),  -omega*sin(nu), 0, sin(nu), cos(nu),  0, omega*cos(nu);
           0,              0,              0, 0,       0,        1, 0]*rel_cylindrical';
    pos = ECEF(1:3);
    velocity = ECEF(4:6);
end

function cost = sum_control_cost(control_traj,dt)
    if iscell(control_traj)
        cost = [multicomplex(0); multicomplex(0); multicomplex(0)];
        for i = 1:length(control_traj)
            cost = cost + dt*[control_traj{i,:}]';
        end
        cost = [repr(cost(1));repr(cost(2));repr(cost(3))];
    else
       cost = [sum(control_traj(:,1));sum(control_traj(:,2));sum(control_traj(:,3))];
    end
end

%{
function kepleriancoord = ECEF_to_keplerian2D(state_vector, mu)
    x = state_vector(1);
    y = state_vector(2);
    z = state_vector(3);
    r = state_vector(1:3);
    v = state_vector(4:6);
    
    h = cross(r,v);
    %p = h^2/mu;
    
    e_vector = cross(v, h)/mu - r/norm(r);
    e = norm(e_vector);
    
    n = cross([0;0;1],h);
    
    i = acos(h(3)/norm(h)); % should = 0
    
    nu = acos(dot(e_vector,r)/(e*norm(r))); %True anomoly
       
    E = 2*atan(tan(nu/2)/sqrt((1+e)/(1-e)));
    
    if n(2)>=0
        argp = acos(n(1)/norm(n));
    else
        argp = 2*pi - acos(n(1)/norm(n));
    end
    if e_vector>=0
        RAAN = acos(dot(n,e_vector)/(e*norm(n)));
    else
        RAAN = 2*pi - acos(dot(n,e_vector)/(e*norm(n)));
    end
    if isnan(RAAN)
        RAAN = 0;
    end
    if  isnan(argp)
        argp = 0;
    end
    
    M = E - e*sin(E);
    
    a = (norm(r)*mu)/(2*mu - norm(v)^2*norm(r));
    %a = h.^2/(mu*(1-e^2));
    kepleriancoord = [a, e, i, argp, RAAN, M];
end


function ECEF2D = keplerian_to_ECEF2D(keplerian_coord, mu, t)
    a = keplerian_coord(1);
    e = keplerian_coord(2);
    i = keplerian_coord(3);
    argp = keplerian_coord(4);
    RAAN = keplerian_coord(5);
    M0 = keplerian_coord(6);
    
    n = sqrt(mu/a^3);
    M = M0 + n*t;
    
    E = keplerEq(M,e);
    
    y = sqrt(1+e)*sin(E/2);
    x = sqrt(1-e)*cos(E/2);
    if x > 0
        nu = -2*atan(y/x)
    elseif x < 0
        if y < 0
            nu = -2*(atan(y/x)-pi);
        else
            nu = -2*(atan(y/x)+pi);
        end
    else
        if y < 0
            nu = -2*(-pi/2);
        elseif y == 0
            nu = 'undefined';
        else
            nu = -2*(pi/2);
        end
    end
    %nu = -2*atan(sqrt((1+e)/(1-e))*tan(E/2));   
    
    r = a*(1-e^2)/(1+e*cos(E));
    
    X = r*(cos(RAAN)*cos(argp+nu)-sin(RAAN)*sin(argp+nu)*cos(i));
    Y = r*(sin(RAAN)*cos(argp+nu)-cos(RAAN)*sin(argp+nu)*cos(i));
    Z = r*(sin(i)*sin(argp+nu));
    
    ECEF2D = [X,Y,Z];
end

function E = keplerEq(M,e)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [rad]. 
    eps = 0.0001;
   	En  = M;
	Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
	while ( abs(Ens-En) > eps )
		En = Ens;
		Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
    end
	E = Ens;
end
%}