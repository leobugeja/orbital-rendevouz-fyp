% Requires Aerospace block set
clc
clear
close all
addpath('eci2orb_gooding','multicomplex')

%% INITIAL VARIABLES
num_orbits = 5;
dt = 10; % in seconds
animation_increments = 10;
reference_altitude = 320e3;
y0 = [40 0 0 0 0 0]; % [x0 y0 z0 xdot0 ydot0 zdot0] [0 0 0 40 0 0]


%RESULTING VARIABLES
mu = 3.986004418e14;
r0 = 6371e3 + reference_altitude;
omega = sqrt(mu/r0^3);
P = (2*pi*1/omega);
t_array = [0:dt:round(P)*num_orbits];
tspan = [0:dt:round(P)*num_orbits];

 
%options = odeset('RelTol',1e-10,'AbsTol', 1e-10);
[t,ode_sol] = ode45(@(t,ode_sol) odefun(t,ode_sol,omega), tspan, y0);

%% Calculate Exact Orbit
%{
ECEF = cylindrical_to_ECEF(y0,t(1),omega,r0); % Earth Centered Earth Fixed
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
exact_sol = zeros(size(t_array,1),6);
earth_spacecraft_orbit = zeros(size(t_array,1),3);
earth_reference_orbit = zeros(size(t_array,1),3);
for k = 1:size(t_array,2)
   exact_sol(k,:) = CW_solution(omega, t_array(k), y0); % relative spacecraft to reference position
   ECEF_spacecraft = cylindrical_to_ECEF(exact_sol(k,:),t_array(k),omega,r0);
   ECEF_reference = cylindrical_to_ECEF([0 0 0 0 0 0],t_array(k),omega,r0);
   earth_spacecraft_orbit(k,:) = ECEF_spacecraft(1:3); % spacecraft orbit - earth centered
   earth_reference_orbit(k,:) = ECEF_reference(1:3); % reference orbit - earth centered
end


%% Euler Solution
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*omega^2 0 0 0 2*omega  0;
     0 0 0 -2*omega 0 0;
     0 0 -omega^2 0 0 0];
I = eye(6);
%abs(eigs(I+A*10))

euler_sol = zeros(length(t_array),6);
euler_sol(1,:) = y0;
euler_sol2 = zeros(length(t_array),6);
euler_sol2(1,:) = y0;
for k = 2:length(t_array)
    euler_sol(k,:) = ( (I + dt*A)*euler_sol(k-1,:)' )';
    euler_sol2(k,:) = ( (I + dt*A)*euler_sol(k-1,:)' )';
end

%% PLOT EARTH SPHERE and Setup Tiles
t = tiledlayout(2,2);
t.TileSpacing = 'none';
t.Padding = 'tight';

%figure('Renderer', 'painters', 'Position', [100 100 1300 900])
hold on
nexttile(1)
%sphere_size = 6371e3;
%[X,Y,Z] = sphere(50);
%surf(X*sphere_size,Y*sphere_size,Z*sphere_size);

%% PLOTTING FIXED LINES
nexttile(1)
hold on
plot3(earth_reference_orbit(:,1),earth_reference_orbit(:,2),earth_reference_orbit(:,3),'r');
plot3(earth_spacecraft_orbit(:,1),earth_spacecraft_orbit(:,2),earth_spacecraft_orbit(:,3),'b');
nexttile(2)
hold on
plot(exact_sol(:,2),exact_sol(:,1),'r');
plot(ode_sol(:,2),ode_sol(:,1),'g');
plot(euler_sol(:,2),euler_sol(:,1),'b');
plot(euler_sol2(:,2),euler_sol2(:,1),'k');

nexttile(4)
hold on
plot(exact_sol(:,4),exact_sol(:,5),'r');
plot(ode_sol(:,4),ode_sol(:,5),'g');
plot(euler_sol(:,4),euler_sol(:,5),'b');
figure(2)
hold on
%plot(earth_reference_orbit(:,2))
%plot(earth_reference_orbit(:,1))
%plot(earth_spacecraft_orbit(:,2))
%plot(earth_reference_orbit(:,2)-earth_spacecraft_orbit(:,2))
plot(100*(exact_sol(:,1)-euler_sol(:,1))./exact_sol(:,1));
plot(100*(exact_sol(:,2)-euler_sol(:,2))./exact_sol(:,2));



%% Plot Settings
figure(2)
legend('euler y error','euler x error')
xlabel('time (s)')
ylabel('% error')
figure(1)
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
ylabel('x_{relative}')
xlabel('y_{relative}')
zlabel('z_{relative}')
x_border = range(exact_sol(:,2))*0.1 + 1;
y_border = range(exact_sol(:,1))*0.1 + 1;
xlim([min(exact_sol(:,2))-x_border max(exact_sol(:,2))+x_border])
ylim([min(exact_sol(:,1))-y_border max(exact_sol(:,1))+y_border])
%legend('Exact HW Solution','ODE Solver HW Solution','Euler Solution','AutoUpdate','off') 
nexttile(3)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
nexttile(4)
%axis equal
xlabel('xdot_{relative}')
ylabel('ydot_{relative}')
zlabel('zdot_{relative}')
x_border = range(exact_sol(:,4))*0.1 + 0.5;
y_border = range(exact_sol(:,5))*0.1 + 0.5;
xlim([min(exact_sol(:,4))-x_border max(exact_sol(:,4))+x_border])
ylim([min(exact_sol(:,5))-y_border max(exact_sol(:,5))+y_border])
legend('Exact HW Solution','ODE Solver HW Solution','Euler Solution','AutoUpdate','off') 

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
animated_relative_pos_euler = animatedline('Marker','x','LineWidth',2,'Color','k', 'MaximumNumPoints',1);
nexttile(4)
animated_relative_speed = animatedline('Marker','x','LineWidth',2,'Color','r', 'MaximumNumPoints',1);

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
    addpoints(animated_relative_pos,exact_sol(k,2),exact_sol(k,1));
    addpoints(animated_relative_pos_euler,euler_sol(k,2),euler_sol(k,1));
    nexttile(4)
    addpoints(animated_relative_speed,exact_sol(k,4),exact_sol(k,5));
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


