function runge_kutta_nth_order()
    % Solve n-th order ODE using Runge-Kutta Method
    % The user inputs the equations as a system of first-order ODEs.
    
    % Get user inputs
    n = input('Enter the order of the differential equation (n): ');
    eqns = input('Enter the system of ODEs as a cell array of function handles:\n');
    % Example: {@(t, Y) Y(2), @(t, Y) -2*Y(1) + sin(t)} for y'' = -2y + sin(t)
    
    t0 = input('Enter the initial time t0: ');
    tf = input('Enter the final time tf: ');
    y0 = input('Enter the initial conditions as a vector [y(0), y''(0), ...]: ');
    h = input('Enter the step size h: ');

    % Ensure correct number of equations and initial conditions
    if numel(eqns) ~= n || numel(y0) ~= n
        error('The number of equations and initial conditions must match the order of the ODE.');
    end

    % Initialization
    N = floor((tf - t0) / h); % Number of steps
    t = t0;                  % Start time
    Y = y0(:);               % Initial state vector (column vector)
    time = zeros(N+1, 1);    % Time vector
    Y_sol = zeros(N+1, n);   % Solution matrix

    % Store initial conditions
    time(1) = t;
    Y_sol(1, :) = Y';

    % Runge-Kutta 4th Order Method
    for i = 1:N
        % Compute RK4 coefficients
        k1 = compute_f(eqns, t, Y);
        k2 = compute_f(eqns, t + h/2, Y + h/2 * k1);
        k3 = compute_f(eqns, t + h/2, Y + h/2 * k2);
        k4 = compute_f(eqns, t + h, Y + h * k3);
        
        % Update solution
        Y = Y + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        % Update time
        t = t + h;
        
        % Store results
        time(i+1) = t;
        Y_sol(i+1, :) = Y';
    end

    % Plot Results
    figure;
    for j = 1:n
        subplot(n, 1, j);
        plot(time, Y_sol(:, j), 'LineWidth', 1.5);
        xlabel('Time t');
        ylabel(['y^{' num2str(j-1) '}']);
        grid on;
        title(['Solution for y^{' num2str(j-1) '}(t)']);
    end
end

function f_vec = compute_f(eqns, t, Y)
    % Helper function to compute the vector of derivatives
    n = numel(eqns);
    f_vec = zeros(n, 1);
    for i = 1:n
        f_vec(i) = eqns{i}(t, Y);
    end
end
