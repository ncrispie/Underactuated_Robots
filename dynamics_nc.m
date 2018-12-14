load('dynamics_params.mat');

N = 101; % number of points, meaning N-1 windows
T = 10;
dt = T / (N-1);

% pick cost constants
% we're leaving our various other costs: 
% q_N, q_bold_N, q_n, q_bold_n, r_n, P_n
% test
% Q_n = diag([1 1 1 1 1 1]);
Q_n = diag([0 0 0 10 10 10]);
Q_n(1:3,1:3) = ones(3,3);
Q_N = 10*Q_n;
q_n = 1;
q_bold_n = [-2 ; -2 ; -2 ; 0 ; 0 ; 0];
R = 0.01*diag([1 1 1]);

% theta = 0 is straight up. theta positive is clockwise.

% get symbolic derivatives
[x, u, x_dot_nonlin_sym] = get_dyn2();
A_sym = jacobian(x_dot_nonlin_sym, x);
B_sym = jacobian(x_dot_nonlin_sym, u);

% lambda accessor functions
x_dot_nonlin = @(x_eval, u_eval) double(subs(x_dot_nonlin_sym,  [x ; u], [x_eval ; u_eval]));
x_dot_lin = @(x_about, u_about, x_eval, u_eval) ...
    double(subs(taylor(x_dot_nonlin_sym, [x ; u], [x_about ; u_about],'Order', 2), ...
    [x ; u], [x_eval ; u_eval]));

% continuous to discrete
% see https://en.wikibooks.org/wiki/Control_Systems/Digital_State_Space#Discrete_Coefficient_Matrices
A_cont = @(x_about, u_about) double(subs(A_sym, [x ; u], [x_about ; u_about]));
B_cont = @(x_about, u_about) double(subs(B_sym, [x ; u], [x_about ; u_about]));
A = @(x_about, u_about) expm(dt*A_cont(x_about, u_about));
B = @(x_about, u_about) A_cont(x_about, u_about) \ (A(x_about, u_about) - eye(6))*B_cont(x_about, u_about);

% usage
% remember, linearize about the nominal trajectory
% x_dot_lin(x_0, u_0, x_0, u_0);
% x_dot_nonlin_0 = x_dot_nonlin(x_0, u_0);
% A_0 = A(x_0, u_0);
% B_0 = B(x_0, u_0);

% generate desired trajectory
x_nominal = zeros(6, N);
u_nominal = zeros(3, N);
x_nominal(:, 1) = [0.01 ; pi+0.01 ; pi/2 + 0.01 ; 0 ; 0 ; 0];
x_desired = [-pi/3 ; pi/2 ; pi/3 ; 0 ; 0 ; 0];
K = [10*eye(3) 2*eye(3)];
for n = 1:N-1
    disp(n)
    
    % get control
    control = K * (x_desired - x_nominal(:, n));
    % control(control > 10) = 10; control(control < -10) = -10;
    u_nominal(:, n) = control;
    
    % 2nd order runge kutta
    % refer to http://lpsa.swarthmore.edu/NumInt/NumIntSecond.html
    k1 = x_dot_nonlin(x_nominal(:, n), u_nominal(:, n));
    x_mid = x_nominal(:, n) + k1*dt/2;
    k2 = x_dot_nonlin(x_mid, u_nominal(:, n));
    x_nominal(:, n+1) = x_nominal(:, n) + dt*k2;
end

%BIG LOOP START

big_loop_its = 10;
x_nominal_history = zeros(6, N, big_loop_its+1);
u_nominal_history = zeros(3, N, big_loop_its+1);
merit_history = zeros(1, big_loop_its+1);
for big_looper = 1:big_loop_its
    % begin forward pass
    A_n = zeros(6,6,N);
    B_n = zeros(6,3,N);
    A_n(:,:,N) = A(x_nominal(:,N), u_nominal(:,N));
    B_n(:,:,N) = B(x_nominal(:,N), u_nominal(:,N));
    for n = N-1:-1:1
        disp(n)
        A_n(:,:,n) = A(x_nominal(:, n), u_nominal(:, n));
        B_n(:,:,n) = B(x_nominal(:, n), u_nominal(:, n));

        % eqn 17
    end

    % big S is S2 (8x8)
    % little s bold is S1 (8x1)
    % little s is S0 (scalar)

    % eqn 29
    S2 = zeros(6,6,N);
    S1 = zeros(6,1,N);
    S0 = zeros(N,1);
    S2(:,:,N) = Q_N;
    S1(:,N) = q_bold_n;
    S0(N) = q_n;
    g_n = zeros(3, 1, N);
    G_n = zeros(3, 6, N);
    H_n = zeros(3, 3, N);
    for n = N-1:-1:1
        disp(n)
        g_n(:,:,n) = B_n(:,:,n)' * S1(:,n+1);
        G_n(:,:,n) = B_n(:,:,n)' * S2(:,:,n+1) * A_n(:,:,n);
        H_n(:,:,n) = R + B_n(:,:,n)' * S2(:,:,n+1) * B_n(:,:,n);

        S2(:,:,n) = Q_n + A_n(:,:,n)'*S2(:,:,n+1)*A_n(:,:,n) - G_n(:,:,n).'*pinv(H_n(:,:,n))*G_n(:,:,n);
        S1(:,n) = q_bold_n + A_n(:,:,n)' * S1(:,n+1) - G_n(:,:,n)'*pinv(H_n(:,:,n))*g_n(:,:,n);
        S0(n) = q_n + S0(n+1) + - (1/2)*g_n(:,:,n).'*pinv(H_n(:,:,n))*g_n(:,:,n);
    end

    %L: [3 6]
    % l [3 1]
    %H [3 3]
    % h [ 3 1]

    %line search
    sigma = 1000;
    u_candidate = zeros(3,N);
    x_candidate = zeros(6,N);
    x_candidate(:,1) = x_nominal(:,1); %preseed with initial value;
    merit0 = 0;

    %first get merit function for alpha = 0
    alpha = 0;
    for n = 1:1:N-1
        disp(n)
        u_candidate(:,n) = u_nominal(:,n) - alpha * (g_n(:,:,n) + G_n(:,:,n) * (x_candidate(:,n) - x_nominal(:,n)));  
        k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
        x_mid = x_candidate(:, n) + k1*dt/2;
        k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
        x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
        merit0 = merit0 + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
            + u_candidate(:,n)'*R*u_candidate(:,n);
    end

    %get merit function for alpha = 1
    alpha = 1;
    merit_alpha = 0;
    for n = 1:1:N-1
        disp(n)
        u_candidate(:,n) = u_nominal(:,n) - alpha * (g_n(:,:,n) + G_n(:,:,n) * (x_candidate(:,n) - x_nominal(:,n)));  
        k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
        x_mid = x_candidate(:, n) + k1*dt/2;
        k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
        x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
        merit_alpha = merit_alpha + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
            + u_candidate(:,n)'*R*u_candidate(:,n);
    end

    %compare and evaluate futher if necessary
    for iters = 1:13
        disp(iters)
        if merit_alpha < merit0
            break %if find a merit less than previous merit, break
        else
            merit_alpha= 0;
            alpha = alpha * 0.7;
            %find a new merit
            for n = 1:1:N-1
                u_candidate(:,n) = u_nominal(:,n) - alpha * (g_n(:,:,n) + G_n(:,:,n) * (x_candidate(:,n) - x_nominal(:,n)));  
                k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
                x_mid = x_candidate(:, n) + k1*dt/2;
                k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
                x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
                merit_alpha = merit_alpha + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
                    + u_candidate(:,n)'*R*u_candidate(:,n);
            end
        end
    end
    
    % add to history
    x_nominal_history(:, :, big_looper) = x_nominal;
    u_nominal_history(:, :, big_looper) = u_nominal;
    merit_history(big_looper) = merit0;
    
    u_nominal = u_candidate;
    x_nominal = x_candidate;
    
    disp(big_looper);
end

x_nominal_history(:, :, end) = x_nominal;
u_nominal_history(:, :, end) = u_nominal;
merit_history(end) = merit_alpha;