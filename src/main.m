clear
clc
close all force
addpath('multicomplex')

back_track_line_method = "newton_raphson"; % or "gradient_descent"
sensitivity_method = @fdTrajectory; % for finite difference method, alternatively use @mxTrajectory for multicomplex step method
dynamics_function = @linearEq; % or @nonlinearEq
objective_function = "radius_objective";
initial_variables =  [0.001, 0.001, -0.001]; % [a_x, a_y, a_z]
initial_state = [4000 0 0 0 -2 2]; % [x y z x_velocity y_velocity z_velocity]
reference_altitude = 400e3; % meters above sea level
integration_time_step = 2; % in seconds
numb_of_integration_steps = 2000;
sensitivty_step_size = 1e-7; % for finite difference and use 1e-15 for multicomplex step
calculate_hessian = true; % false if only want to calculate gradient

[control_parameters,total_fuel] = backTrackLineOptimisation(back_track_line_method, sensitivity_method, dynamics_function, objective_function, initial_variables, initial_state, reference_altitude, integration_time_step, numb_of_integration_steps, sensitivty_step_size, calculate_hessian);

%[control_parameters,total_fuel] = fminconOptimisation(@mxTrajectory,@linearEq, "lagrange_fuel_radius_final_burn",  [0.001, 0.001, -0.001], [4000 0 0 0 -2 2],  400e3, 2, 2000, 1e-15, true)
%[control_parameters,total_fuel] = fminconOptimisation(@fdTrajectory,@linearEq, "lagrange_fuel_radius_final_burn",  [0.001, 0.001, -0.001], [0 10000 0 0 0 0],  400e3, 2, 2000, 1e-7, true)

% Two Relatvie Trajectory Plots:
% plotTrajectories([-0.007698326880490,-6.657344256249828e-04,-9.379846073813817e-04,7.628337504630399e+03/2000], @linearEq, "lagrange_fuel_radius_final_burn", [0.001, 0.001, -0.001], [4000 0 0 0 -2 2], 400e3, 7628/2000, 2000, 1e-7, true, "3d_plots");

% Earth Centered Trajectory Plot:
% plotTrajectories([0.006683775731408,0.001481066172619,-3.068838511224377e-04,2.726028989948144e+03/2000], @linearEq, "lagrange_fuel_radius_final_burn", [0.001, 0.001, -0.001], [-4000 -10000 0 0 3 10], 400e3, 2.726028989948144e+03/2000, 2000, 1e-7, true, "ECEF_position");

% Sensitivity Vector Trajectory Plot:
% plotTrajectories([-0.007698326880490,-6.657344256249828e-04,-9.379846073813817e-04,7.628337504630399e+03/2000], @linearEq, "lagrange_fuel_radius_final_burn", [0.001, 0.001, -0.001], [4000 0 0 0 -2 2], 400e3, 2.726028989948144e+03/2000, 2000, 1e-7, true, "sensitivity_quiver");

% Sensitvity Accuracy Comparison Plot:
% plotSensitivityAccuracies(@nonlinearEq, 'radius_objective', [0.001, 0.001, -0.001], [-4000 -10000 0 0 3 10], 400e3, 5, 1000, true)


function [variables,fval]= fminconOptimisation(sensitivity_method, dynamics_fun, objective_function, initial_variables,t0_state, ref_altitude, intial_time_step, num_steps, h,use_hessian)
% Function which iteratively optimises the trajectory using fmincon and for various different types of objective functions

    clear historical
    lb = [-.1;-.1;-1;0];
    ub = [.1;.1;1;20];
    if objective_function == "lagrange_fuel_radius_final_burn"
        options = optimoptions('fmincon','Display','iter','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn', @(variables, lambda)hessianfcn(variables, lambda, "lagrange_fuel_radius_final_burn")); % ,'SpecifyConstraintGradient',true
        objective = @(variables)sensitivity_method(dynamics_fun, "lagrange_fuel_radius_final_burn", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables , intial_time_step]; %
        [variables,fval] = fmincon(objective,initial_paramaters,[],[],[],[],lb,ub,@(variables)nonlcon(variables,"radius_constraint"),options);
    
    elseif objective_function == "lagrange_fuel_radius"
        options = optimoptions('fmincon','Display','iter','OutputFcn',@outfun,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn', @(variables, lambda)hessianfcn(variables, lambda, "lagrange_fuel_radius")); % ,'SpecifyConstraintGradient',true
        objective = @(variables)sensitivity_method(dynamics_fun, "lagrange_fuel_radius", variables,  t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian);
        initial_paramaters = [initial_variables , intial_time_step]; %
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
        global global_correction_fuel
        stop = false;

        switch state
            case 'init'
                uiFig = uifigure('Position',[100 100 1525 475]);
                uT = uitable(uiFig);
                uT.ColumnName = {'iter', 'fun calls', 'step size', 'objective', 'trajectory fuel', 'correction fuel', 'radius', 'x_f', 'y_f', 'z_f', 'velocity radius', 'x_dot','y_dot','z_dot', 'time','ax', 'ay', 'az'};
                uT.Units='Normalized';
                uT.Position = [.05 .05 .9 .9];
                drawnow
                %plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, state)
                plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, "init_3d_overlay")
            case 'iter'
                if size(x,2) == 4
                    dt = x(end);
                else
                    dt = intial_time_step;
                end
                uT.Data(end+1,:) = [optimValues.iteration optimValues.funccount optimValues.stepsize optimValues.fval global_traj_fuel global_correction_fuel global_radius global_final_state(1:3) global_velocity_radius global_final_state(4:6) dt*num_steps x(1:3)];
                %plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, state)
                plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian, "iter_3d_overlay")
                drawnow
                
                if global_radius < 1 % && global_velocity_radius < 0.1
                    %plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian,state)
                    save('temp','uT')
                    stop = true;
                end
            case 'done'
                %plotTrajectories(x, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian)
                drawnow
            otherwise
        end
    end
end

function [objective, independent_variables, final_state] = backTrackLineOptimisation(method, sensitivity_method, dynamics_fun,objective_function, control_parameters, t0_state, ref_altitude, dt, num_steps, h, use_hessian)
% Back track line method for gradient descent and Newton Raphson methods

    start_totaltime = tic;
    uiFig = uifigure('Position',[100 100 1725 475]);
    uT = uitable(uiFig);
    uT.ColumnName = {'iter', 'fun calls', 'step size', 'objective', 'trajectory fuel', 'correction fuel', 'radius', 'x_f', 'y_f', 'z_f', 'x_dot','y_dot','z_dot', 'time','ax', 'ay', 'az'};
    uT.Units='Normalized';
    uT.Position = [.05 .05 .9 .9];
    drawnow
    
    limit = 1e-1;
    beta = 0.3;
    alpha = 0.1;  
    
    independent_variables = [control_parameters, dt];
    [objective, gradient, hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables);
    
    if size(independent_variables,2) == 4
       dt =  independent_variables(4);
    end
    if method == "gradient_descent"
        hessian = eye(3);
        default_t = 1e-13;
    else
        default_t = 1;
    end
    
    plotTrajectories(independent_variables, dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, dt, num_steps, h, use_hessian, 'init')
    uT.Data(1,:) = [0 1 0 objective fuel sum(abs(final_state(4:6))) radius final_state(1:6) dt*num_steps independent_variables(1:3)];
    drawnow
    f_calls = 1;
    
    
    for i = 1:50
        t = default_t;
        local_calls = 1;
        f_calls = f_calls + 1;
        %Backtracking line search
        while true
            [offset_objective, offset_gradient, offset_hessian, final_state, radius, fuel] = calc_obj_grad(independent_variables-t*gradient'/hessian); % /hessian

            
            if offset_objective > objective - alpha*t*gradient'/hessian*gradient && local_calls < 10 && offset_objective > limit %  alpha*t*gradient'/hessian*gradient && local_calls < 10 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
                local_calls = local_calls + 1;
            else
                break
            end
        end
        
        objective = offset_objective;
        independent_variables = independent_variables-t*gradient'/hessian %/hessian; %2e-16
        gradient = offset_gradient;
        if method == "newton_raphson";
            hessian = offset_hessian;
        end
        if size(independent_variables,2) == 4
            dt =  independent_variables(4);
        end
        plotTrajectories(independent_variables, dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, dt, num_steps, h, use_hessian, 'iter')
        uT.Data(1+i,:) = [i f_calls t objective fuel sum(abs(final_state(4:6))) radius final_state(1:6) dt*num_steps independent_variables(1:3)];
        drawnow
        if or(local_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    end_totaltime = toc(start_totaltime)
    
    function [obj, grad, hessian, state, radius, fuel] = calc_obj_grad(independent_variables)
        [obj, grad, hessian, state, radius, fuel] = sensitivity_method(dynamics_fun, objective_function, independent_variables, t0_state, ref_altitude, dt, num_steps, h, use_hessian);
    end

end

function [objective_value, fun_gradient, hessian, final_state, radius, fuel] = mxTrajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian)
% Function which returns the multicomplex step senstivities for given parameters

    global global_radius
    global global_traj_fuel
    global global_final_state
    global global_velocity_radius
    global global_correction_fuel
    
    persistent historical
    
    if objective_function == "plot_trajectory" || isempty(historical) || h~=historical.control_step || any(historical.variables(1:3) ~= variables(1:3)) || (size(historical.variables,2)==4 && (size(variables,2) == 4 && historical.variables(4) ~= variables(4)))
        
        control_parameters = variables(1:3);
        linear_eq = @linearEq;
        if size(variables,2) == 4
            time_step = variables(end);
            if isequal(dynamics_fun, linear_eq)
                if use_hessian
                    dt = multicomplex(inputconverter(time_step,[5 6],h));
                else
                    dt = multicomplex(inputconverter(time_step,[5],h));
                end
            else
                if use_hessian
                    dt = multicomplex(inputconverter(time_step,[7 8],h));
                else
                    dt = multicomplex(inputconverter(time_step,[4],h));
                end
            end
        else
            dt = time_step;
        end

        historical.variables = variables;
        historical.time_step = time_step;
        historical.control_step = h;

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
            if use_hessian
                control_perturbation.ax = multicomplex(inputconverter(0,[1 3],h));
                control_perturbation.ay = multicomplex(inputconverter(0,[2 4],h));
                control_perturbation.az = multicomplex(inputconverter(0,[1 2],h));
            else
                control_perturbation.ax = multicomplex(inputconverter(0,1,h));
                control_perturbation.ay = multicomplex(inputconverter(0,2,h));
                control_perturbation.az = multicomplex(inputconverter(0,1,h));
            end
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
        
        % Numerical Propagation
        for k = 2:length(t_array)
            k1 = dynamics_fun(complex_state_traj{k-1},[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds  lin:
            k2 = dynamics_fun(complex_state_traj{k-1}+dt*k1/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014seconds lin:
            k3 = dynamics_fun(complex_state_traj{k-1}+dt*k2/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014 seconds lin:
            k4 = dynamics_fun(complex_state_traj{k-1}+dt*k3,[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds lin:
            complex_state_traj{k} = complex_state_traj{k-1} + dt*(k1 + 2*k2 + 2*k3 + k4)/6; % nonlin:0.008 seconds lin:
            real_state_traj(k,:) = real(complex_state_traj{k}); % 0.0001 seconds
        end

        if isequal(dynamics_fun, @linearEq)
            complex_radius.xy = complex_state_traj{end}(1)^2 + complex_state_traj{end}(2)^2;
            complex_radius.z = complex_state_traj{end}(3)^2;
        else
            complex_radius = complex_state_traj{end}(1)^2 + complex_state_traj{end}(2)^2 + complex_state_traj{end}(3)^2;
        end
        
           
        [complex_traj_fuel, complex_tot_fuel] = calComplexFuel(dynamics_fun, complex_state_traj, control_traj,dt,num_steps);
        
              
        %Calculate Objective Functions
        historical.traj_fuel = real(dt)*num_steps*(abs(real(control_traj{1,1}))+abs(real(control_traj{1,2}))+abs(real(control_traj{1,3})));
        
        historical.traj_fuel_gradient = mxGradient(complex_traj_fuel, size(variables,2), h,  dynamics_fun)';
        historical.traj_fuel_hessian = mxHessian(complex_traj_fuel, size(variables,2), h, dynamics_fun);
        
        historical.radius = real(complex_state_traj{end}(1))^2+real(complex_state_traj{end}(2))^2+real(complex_state_traj{end}(3))^2;
        historical.radius_gradient = mxGradient(complex_radius, size(variables,2), h,  dynamics_fun)';
        historical.radius_hessian = mxHessian(complex_radius, size(variables,2), h, dynamics_fun);
        
        historical.vel_radius = real_state_traj(end,4)^2+real_state_traj(end,5)^2+real_state_traj(end,6)^2;

        if isequal(dynamics_fun, @linearEq)
            historical.tot_fuel = real(complex_tot_fuel.xy + complex_tot_fuel.z);
        else
            historical.tot_fuel = real(complex_tot_fuel);
        end
        historical.tot_fuel_gradient = mxGradient(complex_tot_fuel, 4, h,  dynamics_fun)';
        historical.tot_fuel_hessian = mxHessian(complex_tot_fuel, 4, h, dynamics_fun);
        
        global_radius = sqrt(historical.radius);
        global_traj_fuel = historical.traj_fuel;
        global_velocity_radius = sqrt(historical.vel_radius);
        global_final_state = real_state_traj(end,:);
        global_correction_fuel = sum(abs(real_state_traj(end,4:6)));
        
        final_state = real_state_traj(end,:); % redundant
        radius = sqrt(historical.radius); % redundant
        fuel =  historical.traj_fuel; % redundant
    end
    
    if objective_function == "fuel_objective"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        
    elseif objective_function == "radius_objective"
        objective_value = historical.radius;
        fun_gradient = historical.radius_gradient;
        hessian = historical.radius_hessian;
        
    elseif objective_function == "velocity_objective"
        objective_value = historical.vel_radius;
        fun_gradient = historical.vel_radius_gradient;
                
    elseif objective_function == "radius_constraint"
        objective_value = historical.radius - 1^2;
        fun_gradient = historical.radius_gradient;
        
    elseif objective_function == "radius_velocity_constraint"
        objective_value = [historical.radius - 1^2, historical.vel_radius - (1/100)^2 ];
        fun_gradient = [historical.radius_gradient,  historical.vel_radius_gradient];
        
    elseif objective_function == "SQP_objective" 
        objective_value = historical.traj_fuel;
        fun_gradient = historical.radius_gradient; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius_final_burn"
        objective_value = historical.tot_fuel;
        fun_gradient = historical.tot_fuel_gradient;
        hessian.objective = historical.tot_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    end
end
    
function [objective_value, fun_gradient, hessian, final_state, radius, fuel] = fdTrajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, fd_control_step, use_hessian)
% Function which returns the finite difference step senstivities for given parameters

    global global_radius
    global global_traj_fuel
    global global_final_state
    global global_velocity_radius
    global global_correction_fuel
    
    persistent historical
    
    if objective_function == "plot_trajectory" || isempty(historical) || fd_control_step~=historical.control_step || any(historical.variables(1:3) ~= variables(1:3)) || size(historical.variables,2)==4 && (size(variables,2) == 4 && historical.variables(4) ~= variables(4)) || size(historical.variables,2) ~= size(variables,2)
        
        control_parameters = variables(1:3);
        if size(variables,2) == 4
           time_step = variables(end);
        end
        
        historical.variables = variables;
        historical.time_step = time_step;
        historical.control_step = fd_control_step;

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
            fd_x.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,1);
            fd_y.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,2);
            fd_z.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,3);
            fd_radius.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,1)^2 + fd_state_traj.(fd_strings(i))(end,2)^2 + fd_state_traj.(fd_strings(i))(end,3)^2;
            fd_traj_fuel.(fd_strings(i)) = dt(i)*num_steps*(abs(fd_control_traj.(fd_strings(i))(1,1)) + abs(fd_control_traj.(fd_strings(i))(1,2)) + abs(fd_control_traj.(fd_strings(i))(1,3)));
            fd_correction_fuel.(fd_strings(i)) = abs(fd_state_traj.(fd_strings(i))(end,4)) + abs(fd_state_traj.(fd_strings(i))(end,5)) + abs(fd_state_traj.(fd_strings(i))(end,6));
            fd_tot_fuel.(fd_strings(i)) = fd_traj_fuel.(fd_strings(i)) + fd_correction_fuel.(fd_strings(i));
            fd_vel_radius.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,4)^2 + fd_state_traj.(fd_strings(i))(end,5)^2 + fd_state_traj.(fd_strings(i))(end,6)^2;
        end

        %Calculate Objective Functions
        historical.traj_fuel = fd_traj_fuel.none;
        historical.traj_fuel_gradient = fdGradient(fd_traj_fuel, size(variables,2), fd_control_step,  dynamics_fun)';
        historical.traj_fuel_hessian = fdHessian(fd_traj_fuel, size(variables,2), fd_control_step, dynamics_fun);
        
        historical.radius = fd_radius.none;
        historical.radius_gradient = fdGradient(fd_radius, size(variables,2), fd_control_step,  dynamics_fun)';
        historical.radius_hessian = fdHessian(fd_radius, size(variables,2), fd_control_step, dynamics_fun);
        
        historical.vel_radius = fd_vel_radius.none;
        historical.vel_radius_gradient = fdGradient(fd_vel_radius, size(variables,2), fd_control_step,  dynamics_fun)';
        historical.vel_radius_hessian = fdHessian(fd_vel_radius, size(variables,2), fd_control_step, dynamics_fun);
        
        historical.tot_fuel = fd_tot_fuel.none;
        historical.tot_fuel_gradient = fdGradient(fd_tot_fuel, size(variables,2), fd_control_step,  dynamics_fun)';
        historical.tot_fuel_hessian = fdHessian(fd_tot_fuel, size(variables,2), fd_control_step, dynamics_fun);
        
        global_radius = sqrt(fd_radius.none);
        global_traj_fuel = fd_traj_fuel.none;
        global_velocity_radius = sqrt(fd_vel_radius.none);
        global_final_state = fd_state_traj.none(end,:);
        global_correction_fuel = sum(abs(fd_state_traj.none(end,4:6)));
        
        final_state = fd_state_traj.none(end,:); % redundant
        radius = sqrt(fd_radius.none); % redundant
        fuel = fd_traj_fuel.none; % redundant
        
    end
       
    if objective_function == "fuel_objective"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        hessian = historical.traj_fuel_hessian;
        
        final_state = global_final_state;
        radius = global_radius;
        fuel = objective_value;
    
    elseif objective_function == "radius_objective"
        objective_value = historical.radius;
        fun_gradient = historical.radius_gradient;
        hessian = historical.radius_hessian;
        
    elseif objective_function == "velocity_objective"
        objective_value = historical.vel_radius;
        fun_gradient = historical.vel_radius_gradient;
                
    elseif objective_function == "radius_constraint"
        objective_value = historical.radius - 1^2;
        fun_gradient = historical.radius_gradient;
        
    elseif objective_function == "radius_velocity_constraint"
        objective_value = [historical.radius - 1^2, historical.vel_radius - (1/100)^2 ];
        fun_gradient = [historical.radius_gradient,  historical.vel_radius_gradient];
        
    elseif objective_function == "SQP_objective" 
        objective_value = historical.traj_fuel;
        fun_gradient = historical.radius_gradient; % Not sure!!
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius_final_burn"
        objective_value = historical.tot_fuel;
        fun_gradient = historical.tot_fuel_gradient;
        hessian.objective = historical.tot_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
        
    elseif objective_function == "lagrange_fuel_radius"
        objective_value = historical.traj_fuel;
        fun_gradient = historical.traj_fuel_gradient;
        hessian.objective = historical.traj_fuel_hessian;
        hessian.constraints = {historical.radius_hessian};
    end
        
    if objective_function == "plot_trajectory"
        objective_value.x = fd_state_traj.none(:,1);
        objective_value.y = fd_state_traj.none(:,2);
        objective_value.z = fd_state_traj.none(:,3);
        objective_value.x_dot = fd_state_traj.none(:,4);
        objective_value.y_dot = fd_state_traj.none(:,5);
        objective_value.z_dot = fd_state_traj.none(:,6);
        fun_gradient = [fdGradient(fd_x, 3, fd_control_step,  dynamics_fun);fdGradient(fd_y, 3, fd_control_step,  dynamics_fun);fdGradient(fd_z, 3, fd_control_step,  dynamics_fun)];
    end
end 


function output = nonlinearEq(x,control_accel,R,mu,omega)
% The nonlinear dynamical equations
    output = [x(4);
              x(5);
              x(6);
              (R+x(1))*(omega+x(5)/R)^2 - mu*(R+x(1))/((R+x(1))^2+x(3)^2)^(3/2) + control_accel(1);
              -2*x(4)*omega*R/(R+x(1)) - 2*x(4)*x(5)/(R+x(1)) + control_accel(2);
              -mu*x(3)/((R+x(1))^2+x(3)^2)^(3/2) + control_accel(3)];
end

function output = linearEq(x,control_accel,~,~,omega)
% The linear dynamical equations (i.e. Clohessy Wiltshire)
   output = [x(4);
       x(5);
       x(6);
       3*omega^2*x(1) + 2*omega*x(5) + control_accel(1);
       -2*omega*x(4) + control_accel(2);
       -omega^2*x(3) + control_accel(3)];
end

function output = cwTrajectory(t0_state, ref_altitude, time_step, num_steps)
% Clohessy Witshire Exact solution trajectory upto time = time_step*num_steps
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
        CW_sol(k,:) = cwSolution(omega, t_array(k), t0_state);
    end
    output.x = CW_sol(:,1);
    output.y = CW_sol(:,2);
    output.z = CW_sol(:,3);
    output.x_dot = CW_sol(:,4);
    output.y_dot = CW_sol(:,5);
    output.z_dot = CW_sol(:,6);
end

function output = cwSolution(omega,t, X0)
% Clohessy Wiltshire Exact Solution for at a time t
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

function [pos,velocity] = cylindricalToECEF(rel_coord,t,omega,r0)
% Convert Cylindrical coordinates to Earth Centered Earth Facing
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

function plotSensitivityAccuracies(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, use_hessian)
% Function to plot comparison between finite difference and multicomplex step accuracy
    fd_r_ax_sensitivity = [];
    cx_r_ax_sensitivity = [];
    fd_r2_axay_sensitivity = [];
    cx_r2_axay_sensitivity = [];
    
    h = 0.1.^[0:50];
    for i = 1:size(h,2) % 50
        disp(h(i))
        [objective_value, fun_gradient, hessian ] = fdTrajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, h(i), use_hessian);
        fd_r_ax_sensitivity(i) = fun_gradient(1);
        fd_r2_axay_sensitivity(i) = hessian(1,2);
        
        [objective_value, fun_gradient, hessian ] = mxTrajectory(dynamics_fun, objective_function, variables, t0_state, ref_altitude, time_step, num_steps, h(i), use_hessian);
        cx_r_ax_sensitivity(i) = fun_gradient(1);
        cx_r2_axay_sensitivity(i) = hessian(1,2);
    end
    
    plot(h, fd_r_ax_sensitivity)
    plot(h, cx_r_ax_sensitivity)
    symlog(gca, 'y', 1)
    set(gca, 'xdir', 'reverse')
    set(gca, 'XScale', 'log')
    hold on
end

function plotTrajectories(final_variables, dynamics_fun, objective_function, initial_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian, state)
% Function to plot final trajectory in various different ways

    persistent final_z_plot final_xy_plot final_xy_dot_plot iter_num colour_map
    green = [0.4660 0.6740 0.1880];
    blue = [0 0.4470 0.7410];
    red = [0.8500 0.3250 0.0980];
    purple = [0.4940, 0.1840, 0.5560];
    yellow = [0.9290 0.6940 0.1250];
    
    if state == "sensitivity_quiver"
        total_num_steps = num_steps;
        
        quiver_vector = zeros(0,3,3);
        trajectory = zeros(0,3);
        for n = 1:80:total_num_steps
            [step_trajectory, quiver_vector(end+1,:,:)] = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, n, h, use_hessian);
            trajectory(end+1,:) = [step_trajectory.x(end),step_trajectory.y(end),step_trajectory.z(end)];
        end
        %[step_trajectory, quiver_vector(end+1,:,:)] = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        %trajectory(end+1,:) = [step_trajectory.x(end),step_trajectory.y(end),step_trajectory.z(end)];
        
        p1 = plot3(step_trajectory.y/1000,step_trajectory.x/1000, step_trajectory.z/1000,'LineWidth',1.2);
        p1.Color(4) = 0.2;
        hold on
        quiver3(trajectory(:,2)/1000,trajectory(:,1)/1000, trajectory(:,3)/1000,quiver_vector(:,2,1)/1000,quiver_vector(:,1,1)/1000, quiver_vector(:,3,1)/1000,1,'Color', yellow,'LineWidth',1,'MaxHeadSize',0.05);
        quiver3(trajectory(:,2)/1000,trajectory(:,1)/1000, trajectory(:,3)/1000,quiver_vector(:,2,2)/1000,quiver_vector(:,1,2)/1000, quiver_vector(:,3,2)/1000,1.3,'k', 'LineWidth',0.9,'MaxHeadSize',0.05);
        quiver3(trajectory(:,2)/1000,trajectory(:,1)/1000, trajectory(:,3)/1000,quiver_vector(:,2,3)/1000,quiver_vector(:,1,3)/1000, quiver_vector(:,3,3)/1000,0.3,'Color', purple,'LineWidth',0.9,'MaxHeadSize',0.03);
        %patch(x,y,t,'edgealpha',0.5,'facecolor','none','edgecolor','flat')
        alpha(0.3) 
        plot3(trajectory(1,2)/1e3,trajectory(1,1)/1e3,trajectory(1,3)/1e3,'o','LineWidth',1.8,'Color',green,'MarkerSize',8)
        plot3(trajectory(end,2)/1e3,trajectory(end,1)/1e3,trajectory(end,3)/1e3,'x','LineWidth',2.2,'Color',red,'MarkerSize',8)
        
        xlabel('$y_{rel} \; (km)$','Interpreter','latex', 'FontSize',14) %, 'Position',[0 -16 -6],'VerticalAlignment','top','HorizontalAlignment','center')
        ylabel('$x_{rel} \; (km)$','Interpreter','latex', 'FontSize',14) %, 'Position',[17 -6 -6],'VerticalAlignment','top','HorizontalAlignment','center')
        zlabel('$z_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)
        legend('Relative Trajectory', 'Sensitivity Vector a_x','Sensitivity Vector a_y', 'Sensitivity Vector a_z','Initial Position','Target')
        grid on
        view(-103, 37)
            
    elseif state == "ECEF_position"
        mu = 3.986004418e14;
        r0 = 6371e3 + ref_altitude;
        omega = sqrt(mu/r0^3);
        
        ECEF_quiv = zeros(0,3);
        ECEF_accel = zeros(0,3);
        ECEF_pos = zeros(0,3);
        target_pos = zeros(0,3);
        
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        for n = 0:(num_steps-1)
            if n>400
                [ECEF_pos(end+1,:), ~] = cylindricalToECEF([final_traj.x(n+1), final_traj.y(n+1), final_traj.z(n+1),final_traj.x_dot(n+1),final_traj.y_dot(n+1),final_traj.z_dot(n+1)], n*time_step, omega,r0);
                [target_pos(end+1,:),~] = cylindricalToECEF(zeros(1,6), n*time_step, omega,r0);
                if rem(n,50) == 0
                    ECEF_quiv(end+1,:) = ECEF_pos(end,:);
                    nu = n*time_step*omega;
                    ECEF_accel(end+1,:) = [cos(nu),   -sin(nu), 0;
                                     sin(nu),    cos(nu), 0
                                     0      ,    0,       1]*[final_variables(1);final_variables(2)*2;final_variables(3)];
                    
                end
            end
                
        end
        
        
        correction_burn = [cos(nu),   -sin(nu), 0;
                               sin(nu),    cos(nu), 0
                               0      ,    0,       1]*[final_traj.x_dot(end);final_traj.y_dot(end);final_traj.z_dot(end)]/2
        hold on
        grid on
        plot3(target_pos(:,2)/1e6,target_pos(:,1)/1e6, target_pos(:,3)/1e3, 'LineWidth',1.2, 'Color',red);
        plot3(ECEF_pos(:,2)/1e6,ECEF_pos(:,1)/1e6, 1.4*ECEF_pos(:,3)/1e4, 'LineWidth',1.2, 'Color', blue);
        quiver3(ECEF_quiv(:,2)/1e6,ECEF_quiv(:,1)/1e6, 1.4*ECEF_quiv(:,3)/1e4,ECEF_accel(:,2),ECEF_accel(:,1),ECEF_accel(:,3)*10, 0.4,'k','LineWidth',0.8)
        quiver3(ECEF_pos(end,2)/1e6,ECEF_pos(end,1)/1e6,1.4*ECEF_pos(end,3)/1e4,correction_burn(2),correction_burn(1),correction_burn(3),'Color',purple,'LineWidth',1.2)
        plot3(ECEF_pos(1,2)/1e6,ECEF_pos(1,1)/1e6,1.4*ECEF_pos(1,3)/1e4,'o','LineWidth',2,'Color',green,'MarkerSize',7)
        plot3(ECEF_pos(end,2)/1e6,ECEF_pos(end,1)/1e6,1.4*ECEF_pos(end,3)/1e4,'x','LineWidth',2,'Color',red,'MarkerSize',7)
        
        sphere_size = 6371e3;
        [X,Y,Z] = sphere(25);
        earth = surf(X*sphere_size/1e6,Y*sphere_size/1e6,Z*sphere_size/1e6);
        set(earth, 'FaceAlpha', 0.1)
        earth.EdgeColor = 'none';
        set(gca, 'XDir','reverse')  

        axis equal
        view([-56 20])
        xlabel('$y \; (\times 10^{3} \, km)$','Interpreter','latex', 'FontSize',14, 'Position',[0 -16 -6],'VerticalAlignment','top','HorizontalAlignment','center')
        ylabel('$x \; (\times 10^{3} \, km)$','Interpreter','latex', 'FontSize',14, 'Position',[17 -6 -6],'VerticalAlignment','top','HorizontalAlignment','center')
        zlabel('$z \; (km)$','Interpreter','latex', 'FontSize',14)
        legend('Target Orbit', 'Chaser Orbit', 'Thrust Vector','Velocity Correction Burn','Start of Trajectory','End of Trajectory', 'Earth (z axis not to scale)')
        yticklabels({'5','0','-5'})
        xticklabels({'5','0','-5','-10'})
        
    elseif state == 'init'
        figure('Renderer', 'painters', 'Position', [400 10 1300 900])
        hold on
        t = tiledlayout(2,2);
        t.TileSpacing = 'none';
        t.Padding = 'tight';
        initial_traj = fdTrajectory(dynamics_fun, "plot_trajectory", initial_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        CW_traj= cwTrajectory(t0_state, ref_altitude, time_step, num_steps);
        
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
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        nexttile(1)
        delete(final_z_plot)
        final_z_plot = plot([final_traj.z;0],[final_traj.z_dot;0],'k');
        nexttile(2)
        delete(final_xy_plot)
        final_xy_plot = plot(final_traj.y,final_traj.x,'k');
        nexttile(4)
        delete(final_xy_dot_plot)
        final_xy_dot_plot = plot([final_traj.y_dot;0],[final_traj.x_dot;0],'k');
    elseif state == 'init_3d_overlay'
        iter_num = 1;
        colour_map = jet(22);
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        
        figure(1)
        hold on
        grid on
        p1 = plot3(final_traj.y/1000,final_traj.x/1000, final_traj.z/1000,'Color',colour_map(iter_num,:),'LineWidth',1.2);
        p1.Color(4) = 0.5; %(iter_num+10)/37;
        Ylm=ylim;
        Xlm=xlim; 
        Zlm=zlim;
        view([16 60])
        xlabel('$y_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)
        ylabel('$x_{rel} \; (km)$','Interpreter','latex', 'FontSize',14,'Position',[22 -4 -5],'VerticalAlignment','top','HorizontalAlignment','center')
        zlabel('$z_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)

        figure(2)
        hold on
        grid on
        p2 = plot3([final_traj.y_dot;0],[final_traj.x_dot;0], [final_traj.z_dot;0],'Color', colour_map(iter_num,:),'LineWidth',1.2);
        p2.Color(4) = 0.5; %(iter_num+10)/37;
        Ylm=ylim;
        Xlm=xlim; 
        Zlm=zlim;
        xlabel('$\dot y_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14)
        ylabel('$\dot x_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14,'Position',[22 -5 -10],'VerticalAlignment','top','HorizontalAlignment','center')
        zlabel('$\dot z_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14)
        view([16 60])
                
    elseif state == 'iter_3d_overlay'
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        figure(1)
        p1 = plot3(final_traj.y/1000,final_traj.x/1000, final_traj.z/1000,'Color',colour_map(iter_num,:),'LineWidth',1.2);
        p1.Color(4) = 0.5; %(iter_num+10)/37;
        
        figure(2)
        p2 = plot3([final_traj.y_dot;0],[final_traj.x_dot;0], [final_traj.z_dot;0],'Color', colour_map(iter_num,:),'LineWidth',1.2);
        p2.Color(4) = 0.5; %(iter_num+10)/37;
        if iter_num < 22
            iter_num = iter_num + 1;
        end
    elseif state == '3d_plots'      
        initial_traj = fdTrajectory(dynamics_fun, "plot_trajectory", initial_variables, t0_state, ref_altitude, 2, num_steps, h, use_hessian);
        final_traj = fdTrajectory(dynamics_fun, "plot_trajectory", final_variables, t0_state, ref_altitude, time_step, num_steps, h, use_hessian);
        CW_traj= cwTrajectory(t0_state, ref_altitude, 2, num_steps);
        
        figure(1)
        grid on
        hold on
        plot3(CW_traj.y/1000,CW_traj.x/1000, CW_traj.z/1000,'Color',red,'LineWidth',1.2) 
        plot3(initial_traj.y/1000,initial_traj.x/1000,initial_traj.z/1000,'Color', purple,'LineWidth',1.2)
        final_z_plot = plot3(final_traj.y/1000,final_traj.x/1000, final_traj.z/1000,'Color',blue,'LineWidth',1.2);
        plot3(t0_state(2)/1000,t0_state(1)/1000,t0_state(3)/1000,'o','LineWidth',2,'Color',green)
        plot3(0,0,0,'x','LineWidth',2,'Color',red)
        
        xlabel('$y_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)
        ylabel('$x_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)
        zlabel('$z_{rel} \; (km)$','Interpreter','latex', 'FontSize',14)
        legend('Zero Control Trajectory', 'Initial Controlled Trajectory','Constant Acceleration Trajectory','Initial Position','Target','AutoUpdate','off') 
        view([-30 40])
        
        figure(2)
        hold on
        grid on

        plot3(CW_traj.y_dot,CW_traj.x_dot, CW_traj.z_dot,'Color',red,'LineWidth',1.2) 
        plot3(initial_traj.y_dot,initial_traj.x_dot,initial_traj.z_dot,'Color',purple,'LineWidth',1.2)
        final_z_plot = plot3([final_traj.y_dot;0],[final_traj.x_dot;0], [final_traj.z_dot;0],'Color', blue,'LineWidth',1.2);
        plot3([final_traj.y_dot(end);0],[final_traj.x_dot(end);0], [final_traj.z_dot(end);0],'k','LineWidth',1.2);
        plot3(t0_state(5),t0_state(4),t0_state(6),'o','LineWidth',2,'Color',green)
        plot3(0,0,0,'x','LineWidth',2,'Color',red)
        
        Ylm=ylim;
        Xlm=xlim;
        Zlm = zlim;
        xlabel('$\dot y_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14,'Position',[mean(Xlm) min(Ylm)*1.8 min(Zlm)],'VerticalAlignment','top','HorizontalAlignment','center')
        ylabel('$\dot x_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14,'Position',[min(Xlm)*1.4 mean(Ylm)*0.8 min(Zlm)],'VerticalAlignment','top','HorizontalAlignment','center')
        zlabel('$\dot z_{rel} \; (m \, s^{-1})$','Interpreter','latex', 'FontSize',14)
        legend('Zero Control Trajectory', 'Initial Parameters Trajectory','Constant Acceleration Trajectory','Velocity Correction Burn','Initial Position','Target','AutoUpdate','off') 
        view([-50 40])
    end 
end


function gradient = mxGradient(obj_fun,gradient_size, h, eqn_type)
% Calculate the gradient for the multicomplex step method
    if isequal(eqn_type, @linearEq)
        gradient = [CXN(obj_fun.xy,1)/h,CXN(obj_fun.xy,2)/h,CXN(obj_fun.z,1)/h];
        if gradient_size == 4
            gradient = [gradient, CXN(obj_fun.xy,5)/h + CXN(obj_fun.z,5)/h]; % check same with obj_fun.z
        end
    else
        gradient = [CXN(obj_fun,1)/h,CXN(obj_fun,2)/h,CXN(obj_fun,3)/h];
        if gradient_size == 4
            if size(obj_fun.zn,2) == 16
                gradient = [gradient, CXN(obj_fun,4)/h];
            else
                gradient = [gradient, CXN(obj_fun,7)/h];
            end
        end
    end
end

function hessian = mxHessian(obj_fun,hessian_size, h, eqn_type)
% Calculate the hessian matrix for the multicomplex step method
    hessian = zeros(hessian_size);
    if isequal(eqn_type, @linearEq) 
        hessian(1,1) = CXN(obj_fun.xy,[1 3])/h^2;
        hessian(2,2) = CXN(obj_fun.xy,[2 4])/h^2;
        hessian(3,3) = CXN(obj_fun.z,[1 2])/h^2;
        
        hessian(1,2) = CXN(obj_fun.xy,[1 2])/h^2;
        %hessian(1,3) = CXN(obj_fun.xy,[1 3])/h^2;
        %hessian(2,3) = CXN(obj_fun.xy,[2 3])/h^2;
        if hessian_size == 4
            hessian(1,4) = CXN(obj_fun.xy,[1 5])/h^2;
            hessian(2,4) = CXN(obj_fun.xy,[2 5])/h^2;
            hessian(3,4) = CXN(obj_fun.z,[1 5])/h^2;
            hessian(4,4) = CXN(obj_fun.xy,[5 6])/h^2 + CXN(obj_fun.z,[5 6])/h^2;
        end
    else
        hessian(1,1) = CXN(obj_fun,[1 4])/h^2;
        hessian(2,2) = CXN(obj_fun,[2 5])/h^2;
        hessian(3,3) = CXN(obj_fun,[3 6])/h^2;
        
        hessian(1,2) = CXN(obj_fun,[1 2])/h^2;
        hessian(1,3) = CXN(obj_fun,[1 3])/h^2;
        hessian(2,3) = CXN(obj_fun,[2 3])/h^2;
        if hessian_size == 4
            hessian(1,4) = CXN(obj_fun,[1 7])/h^2;
            hessian(2,4) = CXN(obj_fun,[2 7])/h^2;
            hessian(3,4) = CXN(obj_fun,[3 7])/h^2;
            hessian(4,4) = CXN(obj_fun,[7 8])/h^2;
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

function [complex_traj_fuel, complex_tot_fuel] = calComplexFuel(dynamics_fun, complex_state_traj, control_traj, dt, num_steps)
% Function to calculate the multicomplex fuel sum

    if isequal(dynamics_fun, @linearEq)
        complex_traj_fuel.xy = multicomplex(0);
        complex_correction_fuel.xy = multicomplex(0);
        for n = 1:2
            if real(control_traj{1,n})>0
                complex_traj_fuel.xy = complex_traj_fuel.xy + control_traj{1,n};
            else
                complex_traj_fuel.xy = complex_traj_fuel.xy - control_traj{1,n};
            end
            if real(complex_state_traj{end}(n+3))>0 % Check if sign of final velocity
                complex_correction_fuel.xy = complex_correction_fuel.xy + complex_state_traj{end}(n+3);
            else
                complex_correction_fuel.xy = complex_correction_fuel.xy - complex_state_traj{end}(n+3);
            end
        end
        complex_traj_fuel.xy = dt*num_steps*complex_traj_fuel.xy;
        complex_tot_fuel.xy = complex_traj_fuel.xy + complex_correction_fuel.xy;

        if real(control_traj{1,3})>0
            complex_traj_fuel.z = control_traj{1,3};
        else
            complex_traj_fuel.z = - control_traj{1,3};
        end
        if real(complex_state_traj{end}(3+3))>0 % Check if sign of final velocity
            complex_correction_fuel.z = complex_state_traj{end}(3+3);
        else
            complex_correction_fuel.z = - complex_state_traj{end}(3+3);
        end
        complex_traj_fuel.z = dt*num_steps*complex_traj_fuel.z;
        complex_tot_fuel.z = complex_traj_fuel.z + complex_correction_fuel.z;
        
    else
        complex_traj_fuel = multicomplex(0);
        complex_correction_fuel = multicomplex(0);
        for n = 1:3
            if real(control_traj{1,n})>0
                complex_traj_fuel = complex_traj_fuel + control_traj{1,n};
            else
                complex_traj_fuel = complex_traj_fuel - control_traj{1,n};
            end
            if real(complex_state_traj{end}(n+3))>0 % Check if sign of final velocity
                complex_correction_fuel = complex_correction_fuel + complex_state_traj{end}(n+3);
            else
                complex_correction_fuel = complex_correction_fuel - complex_state_traj{end}(n+3);
            end
        end
        complex_traj_fuel = dt*num_steps*complex_traj_fuel;
        complex_tot_fuel = complex_traj_fuel + complex_correction_fuel;
    end
end

function gradient = fdGradient(fd_objective, gradient_size, fd_step, eqn_type)
% Calculates the gradient for finite difference method

    gradient = [(fd_objective.ax - fd_objective.axn)/(2*fd_step) , (fd_objective.ay - fd_objective.ayn)/(2*fd_step), (fd_objective.az - fd_objective.azn)/(2*fd_step)];
    if gradient_size == 4
        gradient(end+1) = (fd_objective.dt - fd_objective.dtn)/(2*fd_step);
    end
end

function hessian = fdHessian(fd_objective, hessian_size, fd_step, eqn_type)
% Calculates the hessian matrix for finite difference method

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