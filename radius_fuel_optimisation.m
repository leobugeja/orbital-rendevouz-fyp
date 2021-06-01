clear
clc
close all force
addpath('multicomplex')
clear historical

global constraint_radius
constraint_radius = 1;  % velcoity 0.01

%[objective_value, fun_gradient, hessian, final_state] = complex_trajectory(@nonlin_eqn,"radius_fuel_objective",[-0.005997448690204, 5.044380000798168e-04, 0.030433310357117], 10*[40 -40 100 1 0.1 1], 408e3, 5, 1000, 1e-15, true)
%[objective_value, fun_gradient, hessian, final_state] = fd_trajectory(@nonlin_eqn,[-0.00018, -0.000002, 0.0012], [40 -400 500 0 0 1], 408e3, 5, 1000, 1e-7)

%[control_parameters,final_radius] = OptimiseTrajectory(@complex_trajectory, @nonlin_eqn,'radius_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], [40 -400 500 0 0 1],  408e3, 5, 1000, 1e-15, true) %last is if use hessian
%[objective_val, control_parameters, final_state] = BacktrackLineMethod(@complex_trajectory,@nonlin_eqn,'radius_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-20, true)
%[objective_val, control_parameters, final_state] = BacktrackLineMethod(@complex_trajectory,@lin_eq,'radius_fuel_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-20, true)

%[control_parameters,final_radius] = OptimiseTrajectory(@fd_trajectory, @nonlin_eqn, 'radius_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], [40 -400 500 0 0 1],  408e3, 5, 1000, 1e-7, true)
%[objective_val, control_parameters, final_state] = BacktrackLineMethod(@fd_trajectory, @nonlin_eqn,'radius_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)
%[objective_val, control_parameters, final_state] = BacktrackLineMethod(@fd_trajectory, @nonlin_eqn,'radius_time_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)

%[objective_val, control_parameters, final_state] = SQPMethod(@fd_trajectory, @nonlin_eqn,'radius_time_objective',  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)
%[control_parameters,fuel] = OptimiseTrajectory(@fd_trajectory, @nonlin_eqn, 'separate_constraint',  [-0.006, 5.044380000798168e-04, 0.030433310357117], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)
%[control_parameters,fuel] = OptimiseTrajectory(@fd_trajectory, @nonlin_eqn, "lagrange_fuel_radius",  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)
%[control_parameters,fuel] = OptimiseTrajectory(@fd_trajectory, @nonlin_eqn, "lagrange_fuel_radius_velocity",  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 2, 2000, 1e-7, true)
[control_parameters,total_fuel] = OptimiseTrajectory(@fd_trajectory,@lin_eq, "lagrange_fuel_radius",  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 4, 1000, 1e-7, true)


function [variables,fval]= OptimiseTrajectory(sensitivity_method, dynamics_fun, objective_function, initial_variables,t0_state, ref_altitude, intial_time_step, num_steps, h,use_hessian)
    lb = [-.1;-.1;-1;0];
    ub = [.1;.1;1;10];
    if objective_function == "lagrange_fuel_radius_velocity"
        options = optimoptions('fmincon','Display','iter','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn', @(variables, lambda)hessianfcn(variables, lambda, "lagrange_fuel_radius_velocity")); % ,'SpecifyConstraintGradient',true
        objective = @(variables)sensitivity_method(dynamics_fun, "fuel_objective", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables,intial_time_step];
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_velocity_constraint"),options);
        
    elseif objective_function == "lagrange_fuel_radius" % objective includes last correction burn
        options = optimoptions('fmincon','Display','iter','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn', @(variables, lambda)hessianfcn(variables, lambda, "lagrange_fuel_radius")); % ,'SpecifyConstraintGradient',true
        objective = @(variables)sensitivity_method(dynamics_fun, "lagrange_fuel_radius", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables,intial_time_step];
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_constraint"),options);
    
    elseif objective_function == "gradient_fuel_radius"
        options = optimoptions('fmincon','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'Display','iter','SpecifyConstraintGradient',true); % ,'SpecifyConstraintGradient',true
        objective = @(variables)sensitivity_method(dynamics_fun, "fuel_objective", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables,intial_time_step];
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_constraint"),options);
        
    elseif objective_function == "default"
        options = optimoptions('fmincon','OutputFcn',@outfun);
        objective = @(variables)sensitivity_method(dynamics_fun, "fuel_objective", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables,intial_time_step];
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_velocity_constraint"),options);
    
    elseif objective_function == "lagrange_velocity_radius"
        options = optimoptions('fmincon','Display','iter','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn', @(variables, lambda)hessianfcn(variables, lambda, "lagrange_velocity_radius"));
        objective = @(variables)sensitivity_method(dynamics_fun, "velocity_objective", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables,intial_time_step];
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_constraint"),options);

    else
        objective = @(variables)sensitivity_method(dynamics_fun, objective_function, variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        if objective_function == "hessian_radius"
            options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'Display','iter','HessianFcn','objective');
        
        elseif objective_function == "gradient_radius"
            options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter');
        end
            
        [variables,fval] = fmincon(objective,initial_variables,[],[],[],[],lb,ub,[],options);
    end
    
    function hessian = hessianfcn(variables,lambda, constraint_type)
        [~, ~, hessian_struct] = sensitivity_method(dynamics_fun, constraint_type, variables',  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        hessian = hessian_struct.objective;
        for i = 1:size(lambda.ineqnonlin,1)
            hessian = hessian + lambda.ineqnonlin(i)*hessian_struct.constraints{i};
        end
    end

    function  [c,ceq,gradc,gradceq] = nonlcon(variables, constraint_type)
        ceq = [];
        gradceq = [];
        [c, gradc] = sensitivity_method(dynamics_fun, constraint_type, variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
    end

    function stop = outfun(x,optimValues,state)
        persistent uT
        global global_traj_fuel
        global global_radius
        global global_final_state
        global global_velocity_radius
        %global constraint_radius
        global global_correction_fuel
        stop = false;

        switch state
            case 'init'
                uiFig = uifigure('Position',[100 100 1525 475]);
                uT = uitable(uiFig);
                uT.ColumnName ={'iter', 'fun calls', 'step size', 'objective', 'trajectory fuel', 'correction fuel', 'radius', 'x_f', 'y_f', 'z_f', 'velocity radius', 'x_dot','y_dot','z_dot', 'time','ax', 'ay', 'az'};
                uT.Units='Normalized';
                uT.Position = [.05 .05 .9 .9];
                drawnow
                plot_trajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, state)
            case 'iter'
                uT.Data(end+1,:) = [optimValues.iteration optimValues.funccount optimValues.stepsize optimValues.fval global_traj_fuel global_correction_fuel global_radius global_final_state(1:3) global_velocity_radius global_final_state(4:6) x(4)*num_steps x(1:3)];
                %plot_trajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, state)
                drawnow
                
                if global_radius < 1 % && global_velocity_radius < 0.1
                    plot_trajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian,state)
                    stop = true;
                end
            case 'done'
                %plot_trajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian)
                drawnow
            otherwise
        end
    end
end

function plot_trajectories(final_variables, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian, state)
    persistent final_z_plot final_xy_plot final_xy_dot_plot
    
    if state == 'init'
        figure('Renderer', 'painters', 'Position', [400 10 1300 900])
        hold on
        t = tiledlayout(2,2);
        t.TileSpacing = 'none';
        t.Padding = 'tight';
        initial_traj = fd_trajectory(dynamics_fun, "plot_trajectory", initial_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        final_traj = fd_trajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        CW_traj= CW_trajectory(t0_state, ref_altitude, time_step, num_steps);
        
        nexttile(1)
        hold on
        plot(CW_traj.z,CW_traj.z_dot,'r') 
        plot(initial_traj.z,initial_traj.z_dot,'b')
        final_z_plot = plot([final_traj.z;0],[final_traj.z_dot;0],'k');
        plot(0,0,'og','LineWidth',2)

        nexttile(2)
        hold on
        plot(CW_traj.y,CW_traj.x,'r');
        plot(initial_traj.y,initial_traj.x,'b');
        final_xy_plot = plot(final_traj.y,final_traj.x,'k');
        plot(0,0,'og','LineWidth',2)

        nexttile(4)
        hold on
        plot(CW_traj.y_dot,CW_traj.x_dot,'r');
        plot(initial_traj.y_dot,initial_traj.x_dot,'b');
        final_xy_dot_plot = plot([final_traj.y_dot;0],[final_traj.x_dot;0],'k');
        plot(0,0,'og','LineWidth',2)

        %{
        if objective_function == "lagrange_fuel_radius"
            nexttile(1)
            plot([final_traj.z(end) 0], [final_traj.z_dot(end) 0], 'k')
            nexttile(2)
            plot([final_traj.y(end) 0], [final_traj.x(end) 0], 'k')
            nexttile(3)

            nexttile(4)
            plot([final_traj.y_dot(end) 0], [final_traj.x_dot(end) 0], 'k')
        end
        %}

        % Plot Settings
        nexttile(1)
        xlabel('$z_{rel}$ (m)','Interpreter','latex')
        ylabel('$\dot z_{rel}$ (m)','Interpreter','latex')
        title('Relative Z Axis State Space')
        legend('Zero Control Trajectory', 'Initial Controlled Trajectory','Final Controlled Trajectory','Target','AutoUpdate','off') 

        nexttile(2)
        set(gca, 'xdir', 'reverse')
        title('Relative Position Trajectory')
        ylabel('$x_{rel}$ (m)','Interpreter','latex')
        xlabel('$y_{rel}$ (m)','Interpreter','latex')
        set(get(gca,'YLabel'),'Rotation',0)

        nexttile(4)
        title('Relative Velocity Trajectory')
        %view([-10 65])
        %axis equal
        ylabel('$\dot x_{rel}$ (m)','Interpreter','latex')
        xlabel('$\dot y_{rel}$ (m)','Interpreter','latex')
        set(get(gca,'YLabel'),'Rotation',0)
        set(gca, 'xdir', 'reverse')
        legend('Zero Control Trajectory', 'Initial Controlled Trajectory','Final Controlled Trajectory','Target','AutoUpdate','off') 
    
    elseif state == 'iter'
        final_traj = fd_trajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        nexttile(1)
        delete(final_z_plot)
        final_z_plot = plot([final_traj.z;0],[final_traj.z_dot;0],'k');
        nexttile(2)
        delete(final_xy_plot)
        final_xy_plot = plot(final_traj.y,final_traj.x,'k');
        nexttile(4)
        delete(final_xy_dot_plot)
        final_xy_dot_plot = plot([final_traj.y_dot;0],[final_traj.x_dot;0],'k');
    end 
end

function [objective_value, fun_gradient, hessian, final_state, radius, fuel] = complex_trajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian)
    global global_radius
    global global_traj_fuel
    global global_final_state
    global global_velocity_radius
    global global_correction_fuel
    
    persistent historical
    
    if objective_function == "plot_trajectory" || isempty(historical) || any(historical.variables ~= variables) || historical.time_step ~= time_step
        
        linear_eq = @lin_eq;
        if size(variables,2) == 4
            time_step = variables(end);
            if isequal(dynamics_fun, linear_eq)
                dt = multicomplex(inputconverter(time_step,[5 6],h));
            else
                dt = multicomplex(inputconverter(time_step,[7 8],h));
            end
        else
            dt = time_step;
        end
        
        historical.variables = variables;
        historical.time_step = time_step;

        mu = 3.986004418e14;
        r0 = 6371e3 + ref_altitude;
        omega = sqrt(mu/r0^3);
        P = (2*pi/omega);
        for k = 1:num_steps
            t_array(k) = k*dt;
        end

        control_traj = cell(length(t_array),3);
        control_traj(:,:) = {0};
        complex_state_traj = cell(length(t_array),1);
        complex_state_traj{1} = [multicomplex(t0_state(1));multicomplex(t0_state(2));multicomplex(t0_state(3));multicomplex(t0_state(4));multicomplex(t0_state(5));multicomplex(t0_state(6))];
        real_state_traj = zeros(length(t_array),6);
        real_state_traj(1,:) = t0_state;

        if isequal(dynamics_fun, linear_eq)
            control_perturbation.ax = multicomplex(inputconverter(0,[1 3],h));
            control_perturbation.ay = multicomplex(inputconverter(0,[2 4],h));
            control_perturbation.az = multicomplex(inputconverter(0,[1 2],h));
        elseif use_hessian
            control_perturbation.ax = multicomplex(inputconverter(0,[1 4],h));
            control_perturbation.ay = multicomplex(inputconverter(0,[2 5],h));
            control_perturbation.az = multicomplex(inputconverter(0,[3 6],h));
        else
            control_perturbation.ax = multicomplex(inputconverter(0,1,h));
            control_perturbation.ay = multicomplex(inputconverter(0,2,h));
            control_perturbation.az = multicomplex(inputconverter(0,3,h));
        end

        for i = 1:round(length(t_array))
            control_traj{i,1} = control_parameters(1) + control_perturbation.ax;
            control_traj{i,2} = control_parameters(2) + control_perturbation.ay;
            control_traj{i,3} = control_parameters(3) + control_perturbation.az;
        end
        
        start_RK = tic;
        % Numerical Propagation
        for k = 2:length(t_array)
            k1 = dynamics_fun(complex_state_traj{k-1},[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds  lin:
            k2 = dynamics_fun(complex_state_traj{k-1}+dt*k1/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014seconds lin:
            k3 = dynamics_fun(complex_state_traj{k-1}+dt*k2/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014 seconds lin:
            k4 = dynamics_fun(complex_state_traj{k-1}+dt*k3,[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds lin:
            complex_state_traj{k} = complex_state_traj{k-1} + dt*(k1 + 2*k2 + 2*k3 + k4)/6; % nonlin:0.008 seconds lin:
            real_state_traj(k,:) = real(complex_state_traj{k}); % 0.0001 seconds
        end
        %end_RK = toc(start_RK)

        real_radius = real(complex_state_traj{end}(1))^2+real(complex_state_traj{end}(2))^2+real(complex_state_traj{end}(3))^2;
        complex_radius = complex_state_traj{end}(1)^2 + complex_state_traj{end}(2)^2 + complex_state_traj{end}(3)^2;
      
        complex_traj_fuel = multicomplex(0);
        for n = 1:3
            if real(control_traj{1,n})>0
                complex_fuel = complex_fuel + control_traj{1,n};
            else
                complex_fuel = complex_fuel - control_traj{1,n};
            end
        end
        complex_traj_fuel = dt*num_steps*complex_traj_fuel;

        if objective_function == "radius_fuel_objective"
            objective_value = real_radius + real_fuel;
            complex_objective = complex_radius + complex_fuel;
            if isequal(dynamics_fun, linear_eq)
                hessian = complex_hessian(complex_objective,4,h,"lin");
                fun_gradient = [CXn(complex_objective,1)/h,CXn(complex_objective,2)/h,CXn(complex_objective,1)/h,CXn(complex_objective,5)/h];
            else
                hessian = complex_hessian(complex_objective,4,h,"non_lin");
                fun_gradient = [CXn(complex_objective,1)/h,CXn(complex_objective,2)/h,CXn(complex_objective,4)/h,CXn(complex_objective,7)/h];
            end
        else
            objective_value = real_radius;
            complex_objective = complex_radius;
            hessian = complex_hessian(complex_objective,3,h,"non_lin");
            fun_gradient = [CXn(complex_objective,1)/h,CXn(complex_objective,2)/h,CXn(complex_objective,3)/h];
        end

        final_state = real_state_traj(end,:);
        radius = sqrt(real_radius);
        
        %Calculate Objective Functions
        historical.traj_fuel = real(dt)*num_steps*(abs(real(control_traj{1,1}))+abs(real(control_traj{1,2}))+abs(real(control_traj{1,3})));
        historical.traj_fuel_gradient = fd_gradient(fd_traj_fuel, 4, fd_control_step,  dynamics_fun)';
        historical.traj_fuel_hessian = fd_hessian(fd_traj_fuel, 4, fd_control_step, dynamics_fun);
        
        historical.radius = fd_radius.none;
        historical.radius_gradient = fd_gradient(fd_radius, 4, fd_control_step,  dynamics_fun)';
        historical.radius_hessian = fd_hessian(fd_radius, 4, fd_control_step, dynamics_fun);
        
        historical.vel_radius = fd_vel_radius.none;
        historical.vel_radius_gradient = fd_gradient(fd_vel_radius, 4, fd_control_step,  dynamics_fun)';
        historical.vel_radius_hessian = fd_hessian(fd_vel_radius, 4, fd_control_step, dynamics_fun);
        
        historical.tot_fuel = fd_tot_fuel.none;
        historical.tot_fuel_gradient = fd_gradient(fd_tot_fuel, 4, fd_control_step,  dynamics_fun)';
        historical.tot_fuel_hessian = fd_hessian(fd_tot_fuel, 4, fd_control_step, dynamics_fun);
        
        global_radius = sqrt(fd_radius.none);
        global_traj_fuel = fd_traj_fuel.none;
        global_velocity_radius = sqrt(fd_vel_radius.none);
        global_final_state = fd_state_traj.none(end,:);
        global_correction_fuel = sum(abs(fd_state_traj.none(end,4:6)));
        
        final_state = fd_state_traj.none(end,:); % redundant
        radius = fd_radius.none; % redundant
        fuel = fd_traj_fuel.none; % redundant
        
    end
    
    global constraint_radius
    if objective_function == "fuel_objective"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        
    elseif objective_function == "velocity_objective"
        objective_value = historical.vel_radius;
        fun_gradient = historical.vel_radius_gradient;
                
    elseif objective_function == "radius_constraint"
        objective_value = historical.radius - constraint_radius^2;
        fun_gradient = historical.radius_gradient;
        
    elseif objective_function == "radius_velocity_constraint"
        objective_value = [historical.radius - constraint_radius^2, historical.vel_radius - (constraint_radius/100)^2 ];
        fun_gradient = [historical.radius_gradient,  historical.vel_radius_gradient];
        
    elseif objective_function == "SQP_objective" 
        objective_value = historical.traj_fuel;
        fun_gradient = historical.radius_gradient; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius" % INCLUDES FINAL BURN
        objective_value = historical.tot_fuel;
        fun_gradient = historical.tot_fuel_gradient;
        hessian.objective = historical.tot_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
       
    elseif objective_function == "lagrange_fuel_radius_velocity"
        objective_value = nan;
        fun_gradient = nan; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian, historical.vel_radius_hessian};
        
    elseif objective_function == "lagrange_velocity_radius"
        objective_value = nan;
        fun_gradient = nan; % Not sure!!
        hessian.objective = historical.vel_radius_hessian;
        hessian.constraints = {historical.radius_hessian};
    end
end
    
function [objective_value, fun_gradient, hessian, final_state, radius, fuel] = fd_trajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, fd_control_step, use_hessian)
  
    global global_radius
    global global_traj_fuel
    global global_final_state
    global global_velocity_radius
    global global_correction_fuel
    
    persistent historical
    
    if objective_function == "plot_trajectory" || isempty(historical) || any(historical.variables ~= variables) || historical.time_step ~= time_step     
        control_parameters = variables(1:3);
        if size(variables,2) == 4
           time_step = variables(end);
        end
        
        historical.variables = variables;
        historical.time_step = time_step;

        mu = 3.986004418e14;
        r0 = 6371e3 + ref_altitude;
        omega = sqrt(mu/r0^3);
        P = (2*pi/omega);
        for k = 1:num_steps
            t_array(k) = k*time_step;
        end

        fd_strings = string({'none','ax','ay','az','axn','ayn','azn','axay','axaz','ayaz','dt','dtax','dtay','dtaz','dtn'});
        dt = time_step*ones(1,15);
        dt = dt + [0 0 0 0 0 0 0 0 0 0 fd_control_step fd_control_step fd_control_step fd_control_step -fd_control_step];

        for i = 1:size(fd_strings,2)
           fd_state_traj.(fd_strings(i)) = zeros(length(t_array),6);
           fd_state_traj.(fd_strings(i))(1,:) = t0_state;
           if i<=10
               fd_control_traj.(fd_strings(i)) = zeros(length(t_array),3);
           end
        end

        for k = 1:round(length(t_array))
            fd_control_traj.none(k,:) = control_parameters;
            fd_control_traj.ax(k,:) = control_parameters + [fd_control_step, 0, 0];
            fd_control_traj.ay(k,:) = control_parameters + [0, fd_control_step, 0];
            fd_control_traj.az(k,:) = control_parameters + [0, 0, fd_control_step];
            fd_control_traj.axn(k,:) = control_parameters - [fd_control_step, 0, 0];
            fd_control_traj.ayn(k,:) = control_parameters - [0, fd_control_step, 0];
            fd_control_traj.azn(k,:) = control_parameters - [0, 0, fd_control_step];
            fd_control_traj.axay(k,:) = control_parameters + [fd_control_step, fd_control_step, 0];
            fd_control_traj.axaz(k,:) = control_parameters + [fd_control_step, 0, fd_control_step];
            fd_control_traj.ayaz(k,:) = control_parameters + [0, fd_control_step, fd_control_step];
            fd_control_traj.dt(k,:) = fd_control_traj.none(k,:);
            fd_control_traj.dtn(k,:) = fd_control_traj.none(k,:);
            fd_control_traj.dtax(k,:) = fd_control_traj.ax(k,:);
            fd_control_traj.dtay(k,:) = fd_control_traj.ay(k,:);
            fd_control_traj.dtaz(k,:) = fd_control_traj.az(k,:);
        end

        % Numerical Propagation
        for k = 2:length(t_array)
            for i = 1:size(fd_strings,2)
                k1 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:),        [fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
                k2 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt(i)*k1/2,[fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
                k3 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt(i)*k2/2,[fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
                k4 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt(i)*k3,  [fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
                fd_state_traj.(fd_strings(i))(k,:) = fd_state_traj.(fd_strings(i))(k-1,:) + dt(i)*(k1 + 2*k2 + 2*k3 + k4)/6;
            end
        end
        
        for i = 1:size(fd_strings,2)
            fd_radius.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,1)^2 + fd_state_traj.(fd_strings(i))(end,2)^2 + fd_state_traj.(fd_strings(i))(end,3)^2;
            fd_traj_fuel.(fd_strings(i)) = dt(i)*num_steps*(abs(fd_control_traj.(fd_strings(i))(1,1)) + abs(fd_control_traj.(fd_strings(i))(1,2)) + abs(fd_control_traj.(fd_strings(i))(1,3)));
            fd_correction_fuel.(fd_strings(i)) = abs(fd_state_traj.(fd_strings(i))(end,4)) + abs(fd_state_traj.(fd_strings(i))(end,5)) + abs(fd_state_traj.(fd_strings(i))(end,6));
            fd_tot_fuel.(fd_strings(i)) = fd_traj_fuel.(fd_strings(i)) + fd_correction_fuel.(fd_strings(i));
            fd_vel_radius.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,4) + fd_state_traj.(fd_strings(i))(end,5)^2 + fd_state_traj.(fd_strings(i))(1,6)^2;
        end

        %Calculate Objective Functions
        historical.traj_fuel = fd_traj_fuel.none;
        historical.traj_fuel_gradient = fd_gradient(fd_traj_fuel, 4, fd_control_step,  dynamics_fun)';
        historical.traj_fuel_hessian = fd_hessian(fd_traj_fuel, 4, fd_control_step, dynamics_fun);
        
        historical.radius = fd_radius.none;
        historical.radius_gradient = fd_gradient(fd_radius, 4, fd_control_step,  dynamics_fun)';
        historical.radius_hessian = fd_hessian(fd_radius, 4, fd_control_step, dynamics_fun);
        
        historical.vel_radius = fd_vel_radius.none;
        historical.vel_radius_gradient = fd_gradient(fd_vel_radius, 4, fd_control_step,  dynamics_fun)';
        historical.vel_radius_hessian = fd_hessian(fd_vel_radius, 4, fd_control_step, dynamics_fun);
        
        historical.tot_fuel = fd_tot_fuel.none;
        historical.tot_fuel_gradient = fd_gradient(fd_tot_fuel, 4, fd_control_step,  dynamics_fun)';
        historical.tot_fuel_hessian = fd_hessian(fd_tot_fuel, 4, fd_control_step, dynamics_fun);
        
        global_radius = sqrt(fd_radius.none);
        global_traj_fuel = fd_traj_fuel.none;
        global_velocity_radius = sqrt(fd_vel_radius.none);
        global_final_state = fd_state_traj.none(end,:);
        global_correction_fuel = sum(abs(fd_state_traj.none(end,4:6)));
        
        final_state = fd_state_traj.none(end,:); % redundant
        radius = fd_radius.none; % redundant
        fuel = fd_traj_fuel.none; % redundant
        
    end
   
    global constraint_radius
    
    if objective_function == "fuel_objective"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        
    elseif objective_function == "velocity_objective"
        objective_value = historical.vel_radius;
        fun_gradient = historical.vel_radius_gradient;
                
    elseif objective_function == "radius_constraint"
        objective_value = historical.radius - constraint_radius^2;
        fun_gradient = historical.radius_gradient;
        
    elseif objective_function == "radius_velocity_constraint"
        objective_value = [historical.radius - constraint_radius^2, historical.vel_radius - (constraint_radius/100)^2 ];
        fun_gradient = [historical.radius_gradient,  historical.vel_radius_gradient];
        
    elseif objective_function == "SQP_objective" 
        objective_value = historical.traj_fuel;
        fun_gradient = historical.radius_gradient; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius" % INCLUDES FINAL BURN
        objective_value = historical.tot_fuel;
        fun_gradient = historical.tot_fuel_gradient;
        hessian.objective = historical.tot_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
       
    elseif objective_function == "lagrange_fuel_radius_velocity"
        objective_value = nan;
        fun_gradient = nan; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian, historical.vel_radius_hessian};
        
    elseif objective_function == "lagrange_velocity_radius"
        objective_value = nan;
        fun_gradient = nan; % Not sure!!
        hessian.objective = historical.vel_radius_hessian;
        hessian.constraints = {historical.radius_hessian};
    end
        
    if objective_function == "plot_trajectory"
        objective_value.x = fd_state_traj.none(:,1);
        objective_value.y = fd_state_traj.none(:,2);
        objective_value.z = fd_state_traj.none(:,3);
        objective_value.x_dot = fd_state_traj.none(:,4);
        objective_value.y_dot = fd_state_traj.none(:,5);
        objective_value.z_dot = fd_state_traj.none(:,6);
    end
end 

function output = nonlin_eqn(x,control_accel,R,mu,omega)
    output = [x(4);
              x(5);
              x(6);
              (R+x(1))*(omega+x(5)/R)^2 - mu*(R+x(1))/((R+x(1))^2+x(3)^2)^(3/2) + control_accel(1);
              -2*x(4)*omega*R/(R+x(1)) - 2*x(4)*x(5)/(R+x(1)) + control_accel(2);
              -mu*x(3)/((R+x(1))^2+x(3)^2)^(3/2) + control_accel(3)];
end

function output = lin_eq(x,control_accel,~,~,omega)
   output = [x(4);
       x(5);
       x(6);
       3*omega^2*x(1) + 2*omega*x(5) + control_accel(1);
       -2*omega*x(4) + control_accel(2);
       -omega^2*x(3) + control_accel(3)];
end

function gradient = complex_gradient(obj_fun,gradient_size, h, eqn_type)


end

function hessian = complex_hessian(obj_fun,hessian_size, h, eqn_type)
    hessian = zeros(hessian_size);
    if eqn_type == "lin_eqn"
        hessian(1,1) = CXn(obj_fun,[1 2])/h^2;
        hessian(2,2) = CXn(obj_fun,[2 5])/h^2;
        hessian(3,3) = CXn(obj_fun,[3 6])/h^2;
        
        hessian(1,2) = CXn(obj_fun,[1 2])/h^2;
        hessian(1,3) = CXn(obj_fun,[1 3])/h^2;
        hessian(2,3) = CXn(obj_fun,[2 3])/h^2;
        if hessian_size == 4
            hessian(1,4) = CXn(obj_fun,[1 7])/h^2;
            hessian(2,4) = CXn(obj_fun,[2 7])/h^2;
            hessian(3,4) = CXn(obj_fun,[3 7])/h^2;
            hessian(4,4) = CXn(obj_fun,[7 8])/h^2;
        end
    else
        hessian(1,1) = CXn(obj_fun,[1 4])/h^2;
        hessian(2,2) = CXn(obj_fun,[2 5])/h^2;
        hessian(3,3) = CXn(obj_fun,[3 6])/h^2;
        
        hessian(1,2) = CXn(obj_fun,[1 2])/h^2;
        hessian(1,3) = CXn(obj_fun,[1 3])/h^2;
        hessian(2,3) = CXn(obj_fun,[2 3])/h^2;
        if hessian_size == 4
            hessian(1,4) = CXn(obj_fun,[1 7])/h^2;
            hessian(2,4) = CXn(obj_fun,[2 7])/h^2;
            hessian(3,4) = CXn(obj_fun,[3 7])/h^2;
            hessian(4,4) = CXn(obj_fun,[7 8])/h^2;
        end
    end
    hessian(2,1) = hessian(1,2);
    hessian(3,1) = hessian(1,3);
    hessian(3,2) = hessian(2,3);
    if hessian_size == 4
        hessian(4,1) = hessian(1,4);
        hessian(4,2) = hessian(2,4);
        hessian(4,3) = hessian(3,4);
    end
end

function gradient = fd_gradient(fd_objective, gradient_size, fd_step, eqn_type)
    gradient = [(fd_objective.ax - fd_objective.axn)/(2*fd_step) , (fd_objective.ay - fd_objective.ayn)/(2*fd_step), (fd_objective.az - fd_objective.azn)/(2*fd_step)];
    if gradient_size == 4
        gradient(end+1) = (fd_objective.dt - fd_objective.dtn)/(2*fd_step);
    end
end

function hessian = fd_hessian(fd_objective, hessian_size, fd_step, eqn_type)
	hessian = zeros(3);
    hessian(1,1) = (fd_objective.ax - 2*fd_objective.none + fd_objective.axn)/fd_step^2;
    hessian(2,2) = (fd_objective.ay - 2*fd_objective.none + fd_objective.ayn)/fd_step^2;
    hessian(3,3) = (fd_objective.az - 2*fd_objective.none + fd_objective.azn)/fd_step^2;

    hessian(1,2) = (fd_objective.axay - fd_objective.ax - fd_objective.ay + fd_objective.none)/fd_step^2;
    hessian(2,1) = hessian(1,2);
    hessian(1,3) = (fd_objective.axaz - fd_objective.ax - fd_objective.az + fd_objective.none)/fd_step^2;
    hessian(3,1) = hessian(1,3);
    hessian(2,3) = (fd_objective.ayaz - fd_objective.ay - fd_objective.az + fd_objective.none)/fd_step^2;
    hessian(3,2) = hessian(2,3);
    if hessian_size == 4
        hessian(1,4) = (fd_objective.dtax - fd_objective.ax - fd_objective.dt + fd_objective.none)/fd_step^2;
        hessian(4,1) = hessian(1,4);
        hessian(2,4) = (fd_objective.dtay - fd_objective.ay - fd_objective.dt + fd_objective.none)/fd_step^2;
        hessian(4,2) = hessian(2,4);
        hessian(3,4) = (fd_objective.dtaz - fd_objective.az - fd_objective.dt + fd_objective.none)/fd_step^2;
        hessian(4,3) = hessian(3,4);
        hessian(4,4) = (fd_objective.dt - 2*fd_objective.none + fd_objective.dtn)/fd_step^2;
    end
end

function [objective, independent_variables, final_state] = BacktrackLineMethod(sensitivity_method, dynamics_fun,objective_function, control_parameters, t0_state, ref_altitude, dt, num_steps, h, use_hessian)
    start_totaltime = tic;
    uiFig = uifigure('Position',[100 100 1525 475]);
    uT = uitable(uiFig);
    uT.ColumnName ={'iter', 'fun calls', 'line step', 'objective val','radius', 'fuel','x_f', 'y_f', 'z_f', 'dt','ax', 'ay', 'az',};
    uT.Units='Normalized';
    uT.Position = [.1 .1 .8 .8];
    drawnow
    
    limit = 1e-1;
    beta = 0.3;
    alpha = 0.1;    
    if use_hessian
        default_t = 1;
    else
        default_t = 1e-13;
    end
    
    % Minimize Position
    objective_function = "radius_objective";
    independent_variables = control_parameters;
    start_iteration = tic;
    [objective, gradient, hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables);
    end_iteration = toc(start_iteration)
    
    uT.Data(1,:) = [0 1 0 objective radius fuel final_state(1:3) dt independent_variables(1:3)];
    drawnow
    
    for i = 1:15
        t = default_t;
        f_calls = 1;             
        % backtracking line search
        while true
            %start_iteration = tic;
            [offset_objective, offset_gradient, offset_hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables-t*gradient/hessian);
            %end_iteration = toc(start_iteration)
            if offset_objective > objective - alpha*t*gradient/hessian*gradient' && f_calls < 10 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        objective = offset_objective;
        independent_variables = independent_variables - t*gradient/hessian; %2e-16
        gradient = offset_gradient;
        hessian = offset_hessian;
        uT.Data(1+i,:) = [i f_calls t objective radius fuel final_state(1:3) dt independent_variables(1:3)];
        drawnow
        if or(f_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    
    %i=0;
    %t = default_t;
    %f_calls = 1;
    % Minimise Position and Fuel
    independent_variables(end+1) = dt;
    objective_function = "radius_fuel_objective";
    [objective, gradient, hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables);
    uT.Data(i+2,:) = [0 f_calls t objective radius fuel final_state(1:3) dt independent_variables(1:3)];
    drawnow
    for j = 1:30
        t = default_t;
        f_calls = 1;             
        % backtracking line search
        while true
            %start_iteration = tic;
            [offset_objective, offset_gradient, offset_hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables-t*gradient/hessian);
            %end_iteration = toc(start_iteration)
            if offset_objective > objective - alpha*t*gradient/hessian*gradient' && f_calls < 50 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        objective = offset_objective;
        independent_variables = independent_variables - t*gradient/hessian; %2e-16
        gradient = offset_gradient;
        hessian = offset_hessian;
        if objective_function == "radius_fuel_objective"
            dt = independent_variables(end);
        end
        uT.Data(i+j+2,:) = [j f_calls t objective radius fuel final_state(1:3) dt independent_variables(1:3)];
        drawnow
        if or(f_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    end_totaltime = toc(start_totaltime)
    
    function [obj, grad, hessian, state, radius, fuel] = calc_obj_grad(independent_variables)
        if objective_function == "radius_fuel_objective"
            [obj, grad, hessian, state, radius, fuel] = sensitivity_method(dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, independent_variables(end), num_steps, h, use_hessian);
        else
            [obj, grad, hessian, state, radius, fuel] = sensitivity_method(dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, dt, num_steps, h, use_hessian);
        end
    end
end

function [objective, independent_variables, final_state] = SQPMethod(sensitivity_method, dynamics_fun,objective_function, control_parameters, t0_state, ref_altitude, dt, num_steps, h, use_hessian)
    start_totaltime = tic;
    uiFig = uifigure('Position',[100 100 1525 475]);
    uT = uitable(uiFig);
    uT.ColumnName ={'iter', 'fun calls', 'line step', 'lagrangian','lambda','radius', 'fuel','x_f', 'y_f', 'z_f', 'dt','ax', 'ay', 'az',};
    uT.Units='Normalized';
    uT.Position = [.1 .1 .8 .8];
    drawnow
    
    limit = 10;
    beta = 0.3;
    alpha = 0.1;  
    default_t = 1;
    
    % Minimize Position
    objective_function = "radius_objective";
    independent_variables = control_parameters;
    [objective, gradient, hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables);
    
    uT.Data(1,:) = [0 1 0 objective 0 radius fuel final_state(1:3) dt independent_variables(1:3)];
    drawnow
    
    for i = 1:15
        t = default_t;
        f_calls = 1;             
        % backtracking line search
        while true
            %start_iteration = tic;
            [offset_objective, offset_gradient, offset_hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables-t*gradient/hessian);
            %end_iteration = toc(start_iteration)
            if offset_objective > objective - alpha*t*gradient/hessian*gradient' && f_calls < 10 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        objective = offset_objective;
        independent_variables = independent_variables - t*gradient/hessian; %2e-16
        gradient = offset_gradient;
        hessian = offset_hessian;
        uT.Data(1+i,:) = [i f_calls t objective 0 radius fuel final_state(1:3) dt independent_variables(1:3)];
        drawnow
        if or(f_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    % SQP Part
    objective_function = "SQP_objective";
    %independent_variables = control_parameters;  
    independent_variables(end+1) = dt;
    
    lambda = 10;
    
    [objective, gradient, hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables);
    lagrangian = fuel - lambda*radius;
    uT.Data(i+1,:) = [0 1 default_t lagrangian lambda radius fuel final_state(1:3) dt independent_variables(1:3)];
    drawnow
    for j = 1:30
        t = default_t;
        f_calls = 1;
        
        %a =  [-gradient.fuel; -radius];
       
        a = [-gradient.fuel+gradient.radius*lambda ; -radius];
        b = [hessian.fuel-hessian.radius*lambda, -gradient.radius;
             -gradient.radius', 0 ];
        p_vector = b\a;
        
        while true
            %start_iteration = tic;
            [offset_objective, offset_gradient, offset_hessian, final_state, radius, offset_fuel] = calc_obj_grad(independent_variables+t*p_vector(1:4)');
            %end_iteration = toc(start_iteration)
            if offset_fuel > fuel && f_calls < 50 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        fuel = offset_fuel;
        independent_variables = independent_variables+t*p_vector(1:4)';
        lambda = lambda + t*p_vector(5);
        
        gradient = offset_gradient;
        hessian = offset_hessian;
        
        lagrangian = fuel - lambda*radius;
        if objective_function == "SQP_objective"
            dt = independent_variables(end);
        end
        uT.Data(i+j+1,:) = [j f_calls t lagrangian lambda radius fuel final_state(1:3) dt independent_variables(1:3)];
        drawnow
        if or(f_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    end_totaltime = toc(start_totaltime)
    
    function [obj, grad, hessian, state, radius, fuel] = calc_obj_grad(independent_variables)
        if objective_function == "SQP_objective"
            [obj, grad, hessian, state, radius, fuel] = sensitivity_method(dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, independent_variables(end), num_steps, h, use_hessian);
        else
            [obj, grad, hessian, state, radius, fuel] = sensitivity_method(dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, dt, num_steps, h, use_hessian);
        end
    end
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
    output.x = CW_sol(:,1);
    output.y = CW_sol(:,2);
    output.z = CW_sol(:,3);
    output.x_dot = CW_sol(:,4);
    output.y_dot = CW_sol(:,5);
    output.z_dot = CW_sol(:,6);
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