% [0.001, 0.001, -0.001], [-4000 -10000 0 0 3 10] && [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1]
%


% complex nonlin no time, 3 accel :  10.2975, 10.3536, 10.2817
% complex nonlin hessian, 3 accel 2 time:  21.5354
% complex nonlin no time, 6 accel :  60.6217, 61.5252, 61.3622
% complex nonlin hessian, 6 accel 2 time:  770.1504 769.2817 765.2124
% complex nonlin no time, 9 accel :  ~6000

% complex lin, 2 accel:  4.3388, 4.4253, 4.3244
% complex lin, 2 accel, 1 time : 5.2659, 5.2006, 5.2976
% complex lin, 2 accel, 2 time : 6.9196, 6.7835, 6.7427, 6.7646
% complex lin, 4 accel:  6.6036, 6.6395
% complex lin, 4 accel,1 time : 9.9970, 9.9184, 9.8939
% complex lin, 4 accel, 1 time (but second), 20.8581, 21.1240
% complex lin, 4 accel, 2 time: 21.9645 22.2114 22.3577 21.8054 21.0794 22.2487

close all force
%clear
clc
%{
complex42 = [21.9645 22.2114 22.3577 21.8054 21.0794 22.2487];
complex62 = [770.1504 769.2817 765.2124];
complex6 = [];
data = {complex6,complex42,complex62};

imag_parts = [6,8];

for i = 1:size(data,2)
    mean_time(i) = mean(data{i});
    neg_error(i) = min(data{i})-mean(data{i});
    pos_error(i) = max(data{i})-mean(data{i});
end

scatter(imag_parts, mean_time)
set(gca, 'YScale', 'log')
hold on
errorbar(imag_parts,mean_time,neg_error,pos_error,'LineStyle','none')

legend
%}
%% Plot time complexity charts
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
plot([0,3,5,6,8,9],[5,10.3,23,60.2,768.2,6000]./5, 'd', 'MarkerSize',8,'LineWidth',1,'Color', red)
set(gca, 'YScale', 'log')
hold on
plot([0:10],0.6*1.095.^[[0:10].^2]+0.4,'--','LineWidth',0.8,'Color', red)

xlabel('Imaginary Parts')
ylabel('Normalised Run Time')
plot([0,3,5,6,8,9],[0.06,0.12,0.22,0.29,0.47,0.65]./0.06, 'rx', 'MarkerSize',9,'LineWidth',1.2,'Color', blue)
plot([0:10],1.3.^[[0:10]],'--','LineWidth',0.8,'Color', blue)

legend('MX Imaginary Parts vs Time','MX Approximation: y = a^{x\^2}', 'FD Equivalent Computation', 'FD Approximation: y = b^{x}')


%% Plot Sensitivity Accuracy
%{
load('sensitivity_accuracy_results')
%fd_r_ax_sensitivity = fd_r_ax_sensitivity;
%fd_r2_axay_sensitivity = fd_r2_axay_sensitivity;
cx_r_ax_sensitivity = real(cx_r_ax_sensitivity);
cx_r2_axay_sensitivity = real(cx_r2_axay_sensitivity);


r_ax_ref = cx_r_ax_sensitivity(end);
r2_axay_ref = cx_r2_axay_sensitivity(end);

plot(h, abs((fd_r_ax_sensitivity-r_ax_ref)./r_ax_ref),'*-')
hold on
plot(h, abs((cx_r_ax_sensitivity-r_ax_ref)./r_ax_ref),'d-')

plot(h, abs((fd_r2_axay_sensitivity-r2_axay_ref)./r2_axay_ref),'o--','LineWidth',0.8)
plot(h, abs((cx_r2_axay_sensitivity-r2_axay_ref)./r2_axay_ref),'s--')

%symlog(gca, 'y', 1)
set(gca, 'xdir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('FD 1st Order Derivate','MX 1st Order Derivate','FD 2nd Order Derivate','MX 2nd Order Derivate')
xlabel('Step Size h')
ylabel('Error Relative to MX Result at h=10^{-50}')
%}
%% Plot Grad and Newton Raphson
%{
fd_grad.objective = [];
fd_grad.radius = [];
fd_grad.fcalls = [];
fd_grad.traj_fuel = [];
fd_grad.correction_fuel = [];

mx_grad.objective = [];
mx_grad.radius = [];
mx_grad.fcalls = [];
mx_grad.traj_fuel = [];
mx_grad.correction_fuel = [];

fd_newton.objective = [];
fd_newton.radius = [];
fd_newton.fcalls = [];
fd_newton.traj_fuel = [];
fd_newton.correction_fuel = [];

mx_newton.objective = [];
mx_newton.radius = [];
mx_newton.fcalls = [];
mx_newton.traj_fuel = [];
mx_newton.correction_fuel = [];

fd_fmincon.objective = [];
fd_fmincon.radius = [];
fd_fmincon.fcalls = [];
fd_fmincon.traj_fuel = [];
fd_fmincon.correction_fuel = [];

mx_fmincon.objective = [];
mx_fmincon.radius = [];
mx_fmincon.fcalls = [];
mx_fmincon.traj_fuel = [];
mx_fmincon.correction_fuel = [];
%}
%{
clear
load('grad_newton_fmincon_results')

green = [0.4660 0.6740 0.1880];

blue = [0 0.4470 0.7410];
lightblue = [0.1010 0.7450 0.8330];
midblue = [0 0.4470 0.7410];

red = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];
brown = [0.6350 0.0780 0.1840];

%{
plot([1 1000],[1 1],'--', 'Color', green, 'LineWidth', 1.2)
hold on
plot(fd_grad.fcalls,fd_grad.radius,'v-', 'Color',brown, 'MarkerSize',4)
plot(mx_grad.fcalls,mx_grad.radius-1000,'.-.', 'Color', midblue, 'MarkerSize',8)

plot(fd_newton.fcalls,fd_newton.radius,'d-','Color', brown)
plot(mx_newton.fcalls(1:24),mx_newton.radius(1:24),'*-.', 'Color', blue)

plot(fd_fmincon.fcalls,fd_fmincon.radius,'s-','Color', red)
plot(mx_fmincon.fcalls,mx_fmincon.radius,'x-.', 'Color', lightblue)
%}
plot([1 1000],[1 1],'r--', 'LineWidth', 1.2)
hold on
plot(fd_grad.fcalls,fd_grad.radius,'v-', 'MarkerSize',4)
plot(mx_grad.fcalls,mx_grad.radius-1000,'.-.', 'MarkerSize',10)

plot(fd_newton.fcalls,fd_newton.radius,'d-')
plot(mx_newton.fcalls(1:24),mx_newton.radius(1:24),'*-.')

plot(fd_fmincon.fcalls,fd_fmincon.radius,'s-')
plot(mx_fmincon.fcalls,mx_fmincon.radius,'x-.')


text(14,0.6,'Newton Raphson')
text(150,0.6,'MATLAB fmincon')
text(20,25000,'Gradient Descent')



set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylabel('Distance from Target (m)')
xlabel('Number of Function Calls')
legend('Feasibility Threshold','FD Gradient Descent','MX Gradient Descent','FD Newton Raphson','MX Newton Raphson','FD MATLAB fmincon','MX MATLAB fmincon')
%}