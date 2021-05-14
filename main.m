% Requires Aerospace block set
clc
clear
close all
addpath('eci2orb_gooding','multicomplex')

%% INITIAL VARIABLES
num_orbits = 1;
dt = 50; % in seconds
animation_increments = 5;
reference_altitude = 320e3;
y0 = [40 -400 100 0 0 0]; % [x0 y0 z0 xdot0 ydot0 zdot0] [0 0 0 40 0 0]


%RESULTING VARIABLES
mu = 3.986004418e14;
r0 = 6371e3 + reference_altitude;
omega = sqrt(mu/r0^3);
P = (2*pi*1/omega);
t_array = [0:dt:round(P)*num_orbits];
tspan = [0:dt:round(P)*num_orbits];

 
%options = odeset('RelTol',1e-10,'AbsTol', 1e-10);
%[t,ode_sol] = ode45(@(t,ode_sol) odefun(t,ode_sol,omega), tspan, y0);

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

%% CW equation solutions at time t 
CW_sol = zeros(length(t_array),6);
earth_spacecraft_orbit = zeros(length(t_array),3);
earth_reference_orbit = zeros(length(t_array),3);
for k = 1:size(t_array,2)
   CW_sol(k,:) = CW_solution(omega, t_array(k), y0); % relative spacecraft to reference position
   ECEF_spacecraft = cylindrical_to_ECEF(CW_sol(k,:),t_array(k),omega,r0);
   ECEF_reference = cylindrical_to_ECEF([0 0 0 0 0 0],t_array(k),omega,r0);
   earth_spacecraft_orbit(k,:) = ECEF_spacecraft(1:3); % spacecraft orbit - earth centered
   earth_reference_orbit(k,:) = ECEF_reference(1:3); % reference orbit - earth centered
end


%% Euler and Runge-Kutta Methods
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
I = eye(6);
%abs(eigs(I+A*10))

euler_sol = zeros(length(t_array),6);
euler_sol(1,:) = y0;
RK_sol = zeros(length(t_array),6);
RK_sol(1,:) = y0;
tay_sol = zeros(length(t_array),6);
tay_sol(1,:) = y0;

for k = 2:length(t_array)
    %Euler
    euler_sol(k,:) = ( (I + dt*A)*euler_sol(k-1,:)' )';
    %Runge-Kutta
    k1 = dt*A*RK_sol(k-1,:)';
    k2 = dt*A*(RK_sol(k-1,:)' + 0.5*k1);
    k3 = dt*A*(RK_sol(k-1,:)' + 0.5*k2);
    k4 = dt*A*(RK_sol(k-1,:)' + k3);
    RK_sol(k,:) = ( RK_sol(k-1,:)' + 1/6*(k1 + 2*k2 + 2*k3 + k4) )';
    % Fourth Order Taylor expansion
    tay_sol(k,:) = ( (I + dt*A + dt^2*A^2/2 + dt^3*A^3/6 + dt^4*A^4/24)*tay_sol(k-1,:)' )';
end

%% Control Solution
state_traj = zeros(length(t_array),6);
state_traj(1,:) = y0;

control_traj = zeros(length(t_array),3);
for i = 1:length(t_array)
    control_traj(i,1) = -0.00001*multicomplex(inputconverter(1,[1],1e-10));
    control_traj(i,2) = -0.00001*multicomplex(inputconverter(1,[2],1e-10));
    control_traj(i,3) = -0.00001*multicomplex(inputconverter(1,[3],1e-10));
end

STM = (I + dt*A + dt^2*A^2/2 + dt^3*A^3/6 + dt^4*A^4/24); % State Transition Matrix
CTM = (dt*eye(6) + dt^2*A/2 + dt^3*A^2/6 + dt^4*A^3/24); % Control Transition Matrix inv(A)*(STM - eye(6))
for k = 2:length(t_array)
    state_traj(k,:) = STM*state_traj(k-1,:)' + CTM*B*control_traj(k,:)';
end

%% PLOT EARTH SPHERE and Setup Tiles
t = tiledlayout(2,3);
t.TileSpacing = 'none';
t.Padding = 'tight';

%figure('Renderer', 'painters', 'Position', [100 100 1300 900])
hold on
nexttile(1)
sphere_size = 6371e3;
[X,Y,Z] = sphere(50);
surf(X*sphere_size,Y*sphere_size,Z*sphere_size);

%% PLOTTING FIXED LINES
nexttile(1)
hold on
plot3(earth_reference_orbit(:,1),earth_reference_orbit(:,2),earth_reference_orbit(:,3),'r');
plot3(earth_spacecraft_orbit(:,1),earth_spacecraft_orbit(:,2),earth_spacecraft_orbit(:,3),'b');
nexttile(2)
hold on
plot3(CW_sol(:,2),CW_sol(:,1),CW_sol(:,3),'r');
%plot(ode_sol(:,2),ode_sol(:,1),'g');
%plot(euler_sol(:,2),euler_sol(:,1),'b');
%plot(RK_sol(:,2),RK_sol(:,1),'k');
plot3(state_traj(:,2),state_traj(:,1),state_traj(:,3),'b');


nexttile(4)
hold on
plot3(CW_sol(:,4),CW_sol(:,5),CW_sol(:,6),'r');
%plot(ode_sol(:,4),ode_sol(:,5),'g');
%plot(euler_sol(:,4),euler_sol(:,5),'b');
%plot(RK_sol(:,4),RK_sol(:,5),'k');
plot3(state_traj(:,4),state_traj(:,5),state_traj(:,6),'b');

%{
nexttile(5)
hold on
plot3(state_traj(:,1),state_traj(:,2),state_traj(:,3),'b');
plot3(CW_sol(:,1),CW_sol(:,2),CW_sol(:,3),'r');


nexttile(6)
hold on
plot3(state_traj(:,4),state_traj(:,5),state_traj(:,6),'b');
plot3(CW_sol(:,4),CW_sol(:,5),CW_sol(:,6),'r');
%}
%{
figure(2)
hold on
plot(100*(CW_sol(:,1)-euler_sol(:,1))./CW_sol(:,1));
plot(100*(CW_sol(:,2)-euler_sol(:,2))./CW_sol(:,2));
plot(100*(RK_sol(:,1)-RK_sol(:,1))./RK_sol(:,1));
plot(100*(RK_sol(:,2)-RK_sol(:,2))./RK_sol(:,2));
figure(2)
legend('euler y error','euler x error','RK y error','RK x error');
xlabel('time (s)')
ylabel('% error')
figure(1)
%}
%{
figure(3)
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


%% Plot Settings
nexttile(1)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
nexttile(2)
%axis equal
set(gca, 'xdir', 'reverse')
title('Relative Position Trajectory')
ylabel('$x_{rel}$ (m)','Interpreter','latex')
xlabel('$y_{rel}$ (m)','Interpreter','latex')
zlabel('$z_{rel}$ (m)','Interpreter','latex')
set(get(gca,'ZLabel'),'Rotation',0)
x_border = range(state_traj(:,2))*0.1 + 1;
y_border = range(state_traj(:,1))*0.1 + 1;
z_border = range(state_traj(:,3))*0.1 + 1;
%xlim([min(state_traj(:,2))-x_border max(state_traj(:,2))+x_border])
%ylim([min(state_traj(:,1))-y_border max(state_traj(:,1))+y_border])
%zlim([min(state_traj(:,3))-z_border max(state_traj(:,3))+z_border])
view([-10 65])
legend('Zero Control Trajectory', 'Controlled Trajectory','AutoUpdate','off') 
nexttile(3)
axis equal
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('$z$ (m)','Interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0)
nexttile(4)
title('Relative Velocity Trajectory')
view([-10 65])
%axis equal
xlabel('$\dot x_{rel}$ (m)','Interpreter','latex')
ylabel('$\dot y_{rel}$ (m)','Interpreter','latex')
zlabel('$\dot z_{rel}$ (m)','Interpreter','latex')
set(get(gca,'ZLabel'),'Rotation',0)
x_border = range(state_traj(:,4))*0.1 + 0.5;
y_border = range(state_traj(:,5))*0.1 + 0.5;
%xlim([min(state_traj(:,4))-x_border max(state_traj(:,4))+x_border])
%ylim([min(state_traj(:,5))-y_border max(CW_sol(:,5))+y_border])
legend('Zero Control Trajectory', 'Controlled Trajectory','AutoUpdate','off') 
%{
nexttile(5)
xlabel('x')
ylabel('y')
zlabel('z')
x_border = range(CW_sol(:,2))*0.1 + 1;
y_border = range(CW_sol(:,1))*0.1 + 1;
xlim([min(CW_sol(:,2))-x_border max(CW_sol(:,2))+x_border])
ylim([min(CW_sol(:,1))-y_border max(CW_sol(:,1))+y_border])
nexttile(6)
xlabel('xdot_{relative}')
ylabel('ydot_{relative}')
zlabel('zdot_{relative}')
x_border = range(CW_sol(:,4))*0.1 + 0.5;
y_border = range(CW_sol(:,5))*0.1 + 0.5;
xlim([min(CW_sol(:,4))-x_border max(CW_sol(:,4))+x_border])
ylim([min(CW_sol(:,5))-y_border max(CW_sol(:,5))+y_border])
%}
%% PLOT ANIMATED LINE
%{
nexttile(1)
animated_earth_reference_orbit = animatedline('LineWidth',2,'Color','r', 'MaximumNumPoints',10);
animated_earth_spacecraft_orbit = animatedline('LineWidth',2,'Color','b', 'MaximumNumPoints',10);
nexttile(2)
animated_relative_orbit = animatedline('LineWidth',1,'Color','g', 'MaximumNumPoints',10);

for k = 1:2:length(t_array)*t_array(2)
    nexttile(1)
    addpoints(animated_earth_reference_orbit,earth_reference_orbit(k,1),earth_reference_orbit(k,2),earth_reference_orbit(k,3));
    addpoints(animated_earth_spacecraft_orbit,earth_spacecraft_orbit(k,1),earth_spacecraft_orbit(k,2),earth_spacecraft_orbit(k,3));
    nexttile(2)
    addpoints(animated_relative_orbit,orbit_sol(k,1),orbit_sol(k,2));
    drawnow
end
%}


%% PLOTTING ANIMATIONS
figure(1)
nexttile(1)
animated_earth_reference_dot = animatedline('Marker','x','MarkerSize',4,'LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_earth_spacecraft_dot = animatedline('Marker','o','LineWidth',2,'Color','b', 'MaximumNumPoints',1);
nexttile(3)
animated_earth_spacecraft_orbit = animatedline('LineWidth',1,'Color','r', 'MaximumNumPoints',5);
animated_earth_reference_orbit = animatedline('LineWidth',1,'Color','b', 'MaximumNumPoints',5);
animated_earth_spacecraft_dot2 = animatedline('Marker','x','MarkerSize',4,'LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_earth_reference_dot2 = animatedline('Marker','o','MarkerSize',4,'LineWidth',2,'Color','b', 'MaximumNumPoints',1);
nexttile(2)
animated_relative_pos = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_relative_pos_traj = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);
nexttile(4)
animated_relative_speed = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);
animated_relative_speed_traj = animatedline('Marker','x','LineWidth',2,'Color','b', 'MaximumNumPoints',1);

for k = 1:animation_increments:length(t_array)
    nexttile(1)
    addpoints(animated_earth_reference_dot,earth_reference_orbit(k,1),earth_reference_orbit(k,2),earth_reference_orbit(k,3));
    addpoints(animated_earth_spacecraft_dot,earth_spacecraft_orbit(k,1),earth_spacecraft_orbit(k,2),earth_spacecraft_orbit(k,3));
    nexttile(3)
    addpoints(animated_earth_reference_orbit,earth_reference_orbit(k,1),earth_reference_orbit(k,2),earth_reference_orbit(k,3));
    addpoints(animated_earth_spacecraft_orbit,earth_spacecraft_orbit(k,1),earth_spacecraft_orbit(k,2),earth_spacecraft_orbit(k,3));
    addpoints(animated_earth_reference_dot2,earth_reference_orbit(k,1),earth_reference_orbit(k,2),earth_reference_orbit(k,3));
    addpoints(animated_earth_spacecraft_dot2,earth_spacecraft_orbit(k,1),earth_spacecraft_orbit(k,2),earth_spacecraft_orbit(k,3));
    nexttile(2)
    addpoints(animated_relative_pos,CW_sol(k,2),CW_sol(k,1),CW_sol(k,3));
    addpoints(animated_relative_pos_traj,state_traj(k,2),state_traj(k,1),state_traj(k,3));
    nexttile(4)
    addpoints(animated_relative_speed,CW_sol(k,4),CW_sol(k,5),CW_sol(k,6));
    addpoints(animated_relative_speed_traj,state_traj(k,4),state_traj(k,5),state_traj(k,6));
    drawnow
end


t_array = t_array';

%% Functions
function dydt = odefun(t,y,omega)
    dydt = zeros(6,1);
    dydt(1) = y(4);
    dydt(2) = y(5);
    dydt(3) = y(6);
    dydt(4) = 3*omega^2*y(1) + 2*omega*y(5);
    dydt(5) = -2*omega*y(5);
    dydt(6) = -omega^2*y(3);
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

function ECEF = cylindrical_to_ECEF(rel_coord,t,omega,r0)
    nu = t*omega;
    rel_cylindrical = [rel_coord, r0]; % [x y z x_dot y_dot z_dot r0]
    % x_earth = cos(nu)*(r0 + rel_coord(1)) - rel_coord(2)*sin(nu);
    % y_earth = sin(nu)*(r0 + rel_coord(1)) + rel_coord(2)*cos(nu);
    % z_earth = rel_coord(3);
    % x_earth_dot = cos(nu)*(r0 + rel_coord(4)) - rel_coord(5)*sin(nu);
    % y_earth_dot = cos(nu)*(r0 + rel_coord(4)) - rel_coord(5)*sin(nu);
    % z_earth_dot = rel_coord(6);
      
    ECEF = [cos(nu),       -sin(nu),       0, 0,       0,        0, cos(nu);
           sin(nu),        cos(nu),        0, 0,       0,        0, sin(nu);
           0,              0,              1, 0,       0,        0, 0;
           -omega*sin(nu), -omega*cos(nu), 0, cos(nu), -sin(nu), 0, -omega*sin(nu);
           omega*cos(nu),  -omega*sin(nu), 0, sin(nu), cos(nu),  0, omega*cos(nu);
           0,              0,              0, 0,       0,        1, 0]*rel_cylindrical';
end

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

