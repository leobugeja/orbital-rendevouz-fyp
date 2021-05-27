clear
clc
close all force
addpath('multicomplex')

%[objective_value, fun_gradient, hessian, final_state] = complex_trajectory(@nonlin_eqn,[-0.00018, -0.000002, 0.0012], [40 -400 500 0 0 1], 408e3, 5, 1000, 1e-15, false)
%[objective_value, fun_gradient, hessian, final_state] = fd_trajectory(@nonlin_eqn,[-0.00018, -0.000002, 0.0012], [40 -400 500 0 0 1], 408e3, 5, 1000, 1e-7)

%[control_parameters,final_radius] = OptimiseTrajectory(@complex_trajectory, @nonlin_eqn, [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-15, true) %last is if use hessian
[objective_val, control_parameters, final_state] = BacktrackLineMethod(@complex_trajectory,@lin_eq,  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-20, true)

%[control_parameters,final_radius] = OptimiseTrajectory(@fd_trajectory, @nonlin_eqn,  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, false)
%[objective_val, control_parameters, final_state] = BacktrackLineMethod(@fd_trajectory, @nonlin_eqn,  [-9.9997e-4, -6.4424e-4, 0.0078], 10*[40 -40 100 1 0.1 1],  408e3, 5, 1000, 1e-7, true)

function [control_parameters,fval,exitflag,output,lambda,grad,hessian]= OptimiseTrajectory(sensitivity_method, dynamics_fun, initial_control_parameters,t0_state, ref_altitude, intial_time_step, num_steps, h,use_hessian)
    if use_hessian
        options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'Display','iter','HessianFcn','objective');
    else
        options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'Display','iter');
    end
    [control_parameters,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(control_parameters)sensitivity_method(dynamics_fun,control_parameters, t0_state, ref_altitude, intial_time_step, num_steps, h, use_hessian),initial_control_parameters,[],[],[],[],[-.001;-.001;-.01],[.001;.001;.01],[],options);
end

function [objective_value, fun_gradient, hessian, final_state] = complex_trajectory(dynamics_fun, control_parameters, t0_state, ref_altitude, time_step, num_steps, h, use_hessian)
    dt = time_step; %multicomplex(inputconverter(time_step,4,h));
    
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
    
    if use_hessian
        control_perturbation.ax = multicomplex(inputconverter(0,[1 4],h));
        control_perturbation.ay = multicomplex(inputconverter(0,[2 5],h));
        control_perturbation.az = multicomplex(inputconverter(0,[3 6],h));
    else
        control_perturbation.ax = multicomplex(inputconverter(0, [1],h));
        control_perturbation.ay = multicomplex(inputconverter(0, [2],h));
        control_perturbation.az = multicomplex(inputconverter(0, [3],h));
    end
    
    for i = 1:round(length(t_array))
        control_traj{i,1} = control_parameters(1) + control_perturbation.ax;
        control_traj{i,2} = control_parameters(2) + control_perturbation.ay;
        control_traj{i,3} = control_parameters(3) + control_perturbation.az;
    end
    start_RK = tic;
    % Numerical Propagation
    for k = 2:length(t_array)
        clc
        tic
        k1 = dynamics_fun(complex_state_traj{k-1},[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds  lin:0.002
        toc
        tic
        k2 = dynamics_fun(complex_state_traj{k-1}+dt*k1/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014seconds lin:0.006
        toc
        tic
        k3 = dynamics_fun(complex_state_traj{k-1}+dt*k2/2,[control_traj{k,:}],r0,mu,omega); % nonlin:0.014 seconds lin:0.006
        toc
        tic
        k4 = dynamics_fun(complex_state_traj{k-1}+dt*k3,[control_traj{k,:}],r0,mu,omega); % nonlin:0.013 seconds lin:0.0035 
        toc
        tic
        complex_state_traj{k} = complex_state_traj{k-1};% + dt*(k1 + 2*k2 + 2*k3 + k4)/6; % nonlin:0.008 seconds lin: 0.0084
        toc

        real_state_traj(k,:) = real(complex_state_traj{k}); % 0.0001 seconds
    end
    %end_RK = toc(start_RK)
    objective_value = real(complex_state_traj{end}(1))^2+real(complex_state_traj{end}(2))^2+real(complex_state_traj{end}(3))^2;
    
    complex_radius = complex_state_traj{end}(1)^2 + complex_state_traj{end}(2)^2 + complex_state_traj{end}(3)^2;
    fun_gradient = [CXn(complex_radius,1)/CXn(control_perturbation.ax,1),CXn(complex_radius,2)/CXn(control_perturbation.ay,2),CXn(complex_radius,3)/CXn(control_perturbation.az,3)];
    final_state = real_state_traj(end,:);
    
    if use_hessian
        hessian = complex_hessian(complex_radius,3,h);
    else
        hessian = eye(3);
    end
    %{
    for i = 1:size(complex_state_traj{end},1)
        final_repr(i,1) = repr(complex_state_traj{end}(i));
    end
    %}   
end

function [objective_value, fun_gradient, hessian, final_state] = fd_trajectory(dynamics_fun, control_parameters, t0_state, ref_altitude, time_step, num_steps, fd_control_step, use_hessian)
    dt = time_step;
    
    mu = 3.986004418e14;
    r0 = 6371e3 + ref_altitude;
    omega = sqrt(mu/r0^3);
    P = (2*pi/omega);
    for k = 1:num_steps
        t_array(k) = k*dt;
    end
    
    fd_strings = string({'none','ax','ay','az','axn','ayn','azn','axay','axaz','ayaz','dt','dtax','dtay','dtaz'});
    %fd_state_traj = struct;
    %fd_control_traj = struct;
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
        fd_control_traj.dtax(k,:) = fd_control_traj.ax(k,:);
        fd_control_traj.dtay(k,:) = fd_control_traj.ay(k,:);
        fd_control_traj.dtaz(k,:) = fd_control_traj.az(k,:);
        fd_control_traj.dt(k,:) = fd_control_traj.none(k,:);
    end
    
    % Numerical Propagation
    for k = 2:length(t_array)
        for i = 1:size(fd_strings,2)
            k1 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:),        [fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
            k2 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt*k1/2,[fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
            k3 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt*k2/2,[fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
            k4 = dynamics_fun(fd_state_traj.(fd_strings(i))(k-1,:)+dt*k3,  [fd_control_traj.(fd_strings(i))(k-1,:)],r0,mu,omega)';
            fd_state_traj.(fd_strings(i))(k,:) = fd_state_traj.(fd_strings(i))(k-1,:) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
    end
    
    %Calculate Radius's for Objective
    for i = 1:size(fd_strings,2)
       fd_radius.(fd_strings(i)) = fd_state_traj.(fd_strings(i))(end,1)^2 + fd_state_traj.(fd_strings(i))(end,2)^2 + fd_state_traj.(fd_strings(i))(end,3)^2;
    end
    
    %fun_gradient = [(fd_radius.ax - fd_radius.none)/fd_control_step , (fd_radius.ay - fd_radius.none)/fd_control_step, (fd_radius.az - fd_radius.none)/fd_control_step];
    fun_gradient = [(fd_radius.ax - fd_radius.axn)/(2*fd_control_step) , (fd_radius.ay - fd_radius.ayn)/(2*fd_control_step), (fd_radius.az - fd_radius.none)/(2*fd_control_step)];
    if use_hessian
        hessian = zeros(3);
        hessian(1,1) = (fd_radius.ax - 2*fd_radius.none + fd_radius.axn)/fd_control_step^2;
        hessian(2,2) = (fd_radius.ay - 2*fd_radius.none + fd_radius.ayn)/fd_control_step^2;
        hessian(3,3) = (fd_radius.az - 2*fd_radius.none + fd_radius.azn)/fd_control_step^2;

        hessian(1,2) = (fd_radius.axay - fd_radius.ax - fd_radius.ay + fd_radius.none)/fd_control_step^2;
        hessian(2,1) = hessian(1,2);
        hessian(1,3) = (fd_radius.axaz - fd_radius.ax - fd_radius.az + fd_radius.none)/fd_control_step^2;
        hessian(3,1) = hessian(1,3);
        hessian(2,3) = (fd_radius.ayaz - fd_radius.ay - fd_radius.az + fd_radius.none)/fd_control_step^2;
        hessian(3,2) = hessian(2,3);
    else
        hessian = eye(3);
    end
    
    objective_value = fd_radius.none;
    final_state = fd_state_traj.none(end,:);  
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

function hessian = complex_hessian(obj_fun,hessian_size, h)
    %hessian_size = log2(size(obj_fun.zn,2));
    hessian = zeros(hessian_size);
    
    hessian(1,1) = CXn(obj_fun,[1 4])/h^2;
    hessian(2,2) = CXn(obj_fun,[2 5])/h^2;
    hessian(3,3) = CXn(obj_fun,[3 6])/h^2;
    
    hessian(1,2) = CXn(obj_fun,[1 2])/h^2;
    hessian(2,1) = hessian(1,2);
    hessian(1,3) = CXn(obj_fun,[1 3])/h^2;
    hessian(3,1) = hessian(1,3);
    hessian(2,3) = CXn(obj_fun,[2 3])/h^2;
    hessian(3,2) = hessian(2,3);
end

function [objective, control_parameters, final_state] = BacktrackLineMethod(sensitivity_method, dynamics_fun, control_parameters, t0_state, ref_altitude, initial_time_step, num_steps, h, use_hessian)
    start_totaltime = tic;
    uiFig = uifigure('Position',[100 100 1125 475]);
    uT = uitable(uiFig);
    uT.ColumnName ={'iter', 'fun calls', 'step size t', 'objective val','x_f', 'y_f', 'z_f', 'ax', 'ay', 'az',};
    uT.Units='Normalized';
    uT.Position = [.1 .1 .8 .8];
    drawnow
    
    if use_hessian
        default_t = 1;
    else
        default_t = 1e-13;
    end
    
    limit = 1e-10;
    beta = 0.3;
    alpha = 0.1;
    start_iteration = tic;
    [objective, gradient, hessian, final_state] = calc_obj_grad(control_parameters);
    %end_iteration = toc(start_iteration)
    
    uT.Data(1,:) = [0 1 0 objective final_state(1:3) control_parameters(1:3)];
    drawnow
    for i = 1:15
        t = default_t;
        f_calls = 1;             
        % backtracking line search
        while true
            start_iteration = tic;
            [offset_objective, offset_gradient, offset_hessian, final_state] = calc_obj_grad(control_parameters-t*gradient/hessian);
            %end_iteration = toc(start_iteration)
            if offset_objective > objective - alpha*t*gradient/hessian*gradient' && f_calls < 10 && offset_objective > limit
                t = t*beta;
                f_calls = f_calls + 1;
            else
                break
            end
        end
           
        objective = offset_objective;
        control_parameters = control_parameters - t*gradient/hessian; %2e-16
        gradient = offset_gradient;
        hessian = offset_hessian;
        uT.Data(1+i,:) = [i f_calls t objective final_state(1:3) control_parameters(1:3)];
        drawnow
        if or(f_calls >= 10, offset_objective <= limit)
            break
        end
    end
    
    end_totaltime = toc(start_totaltime)
    
    function [obj, grad, hessian, state] = calc_obj_grad(control_parameters)
        [obj, grad, hessian, state] = sensitivity_method(dynamics_fun, control_parameters, t0_state, ref_altitude, initial_time_step, num_steps, h, use_hessian);
    end

end

