load('dynamics_params.mat');

N = 101; % number of points, meaning N-1 windows
T = 10;
dt = T / (N-1);

% pick cost constants
% we're leaving our various other costs: 
% q_N, q_bold_N, q_n, q_bold_n, r_n, P_n
% test
% Q_n = diag([1 1 1 1 1 1]);
% Q_N = 10*Q_n;
Q_n = diag([0 0 0 10 10 10]);
Q_N = 10*Q_n;
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
x_nominal(:, 1) = [pi+0.01 ; -0.01 ; pi+0.01 ; 0 ; 0 ; 0];
x_desired = [2*pi/3 ; -pi/2 ; 2*pi/3 ; 0 ; 0 ; 0];
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

% nonlinear constraint
% g2(x_n,n) = 0
% this also applies at N! g3(x_N,N) = g2(x_n,n)
% ee is end effector
% ee_start = [1 ; -1];
% ee_end = [1 ; 1];
g2 = @(x_n, n) -1 - (L1*cos(x_n(1)) + L2*cos(x_n(2)) + L3*cos(x_n(3)));
g3 = @(x_n, n) [-1 - (L1*cos(x_n(1)) + L2*cos(x_n(2)) + L3*cos(x_n(3))); ...
    1 - (L1*sin(x_n(1)) + L2*sin(x_n(2)) + L3*sin(x_n(3)))];

% linearized constraint
% remember to linearize about the nominal trajectory!
% the x_n passed in here should come from the nominal trajectory
C_n = @(x_n) [L1*sin(x_n(1)) L2*sin(x_n(2)) L3*sin(x_n(3)) 0 0 0];
d_n = @(x_n) -g2(x_n,0);
C_N = @(x_n) [L1*sin(x_n(1))  L2*sin(x_n(2))  L3*sin(x_n(3)) 0 0 0 ; ...
             -L1*cos(x_n(1)) -L2*cos(x_n(2)) -L3*cos(x_n(3)) 0 0 0];
d_N = @(x_n) -g3(x_n,0);

%BIG LOOP START

big_loop_its = 10;
x_nominal_history = zeros(6, N, big_loop_its+1);
u_nominal_history = zeros(3, N, big_loop_its+1);
merit_history = zeros(1, big_loop_its+1);
for big_looper = 1:big_loop_its
    % begin forward pass
    A_n = zeros(6,6,N);
    B_n = zeros(6,3,N);
    M_n = zeros(1,3,N);
    N_n = zeros(1,6,N);
    Proj_n = zeros(3,3,N);
    epsilon_n = zeros(3,1,N);
    U_n = zeros(3,6,N);
    A_tilde_n = zeros(6,6,N);
    B_tilde_n = zeros(6,3,N);
    g_tilde_n = zeros(6,1,N);
    q_tilde_n = zeros(1,1,N);
    q_bold_tilde_n = zeros(6,1,N);
    Q_tilde_n = zeros(6,6,N);
    r_tilde_n = zeros(3,1,N);
    R_tilde_n = zeros(3,3,N);
    P_tilde_n = zeros(3,6,N);
    A_n(:,:,N) = A(x_nominal(:,N), u_nominal(:,N));
    B_n(:,:,N) = B(x_nominal(:,N), u_nominal(:,N));
    for n = N-1:-1:1
        disp(n)
        A_n(:,:,n) = A(x_nominal(:, n), u_nominal(:, n));
        B_n(:,:,n) = B(x_nominal(:, n), u_nominal(:, n));

        M_n(:,:,n) = C_n(x_nominal(:, n+1)) * A_n(:,:,n+1) * B_n(:,:,n);
        N_n(:,:,n) = C_n(x_nominal(:, n+1)) * A_n(:,:,n+1) * A_n(:,:,n);

        % eqn 13
        w = null(M_n(:,:,n));
        Proj_n(:,:,n) = w*inv(w'*w)*w';

        % eqn 15
        epsilon_n(:,:,n) = pinv(M_n(:,:,n))*d_n(x_nominal(:,n));
        U_n(:,:,n) = -pinv(M_n(:,:,n))*N_n(:,:,n);

        % eqn 17
        A_tilde_n(:,:,n) = A_n(:,:,n) + B_n(:,:,n)*U_n(:,:,n);
        B_tilde_n(:,:,n) = B_n(:,:,n) * Proj_n(:,:,n);
        g_tilde_n(:,:,n) = B_n(:,:,n) * epsilon_n(:,:,n);

        q_tilde_n(:,:,n) = epsilon_n(:,:,n)' * R * epsilon_n(:,:,n)/2; % 20
        q_bold_tilde_n(:,:,n) = U_n(:,:,n)' * R * epsilon_n(:,:,n); % 21
        Q_tilde_n(:,:,n) = Q_n + U_n(:,:,n)' * R * U_n(:,:,n); % 22
        r_tilde_n(:,:,n) = Proj_n(:,:,n)*R*epsilon_n(:,:,n); % 23
        R_tilde_n(:,:,n) = Proj_n(:,:,n)*R*Proj_n(:,:,n); % 24
        P_tilde_n(:,:,n) = Proj_n(:,:,n)*R*U_n(:,:,n); % 25
    end

    % big S is S2 (8x8)
    % little s bold is S1 (8x1)
    % little s is S0 (scalar)

    % eqn 29
    S2 = zeros(6,6,N);
    S1 = zeros(6,1,N);
    S0 = zeros(N,1);
    S2(:,:,N) = Q_N;
    S1(:,N) = zeros(6,1);
    S0(N) = 0;
    L_n = zeros(3, 6, N);
    l_n = zeros(3, 1, N);
    for n = N-1:-1:1
        disp(n)
        h = r_tilde_n(:,:,n) + B_tilde_n(:,:,n)' * (S1(:,n+1)+ S2(:,:,n+1)*g_tilde_n(:,:,n));
        G = P_tilde_n(:,:,n) + B_tilde_n(:,:,n)' * S2(:,:,n+1) * A_tilde_n(:,:,n);
        H = R_tilde_n(:,:,n) + B_tilde_n(:,:,n)' * S2(:,:,n+1) * B_tilde_n(:,:,n);
        L_n(:,:,n) = -pinv(H)*G;
        l_n(:,:,n) = -pinv(H)*h;

        S2(:,:,n) = Q_tilde_n(:,:,n) + A_tilde_n(:,:,n)'*S2(:,:,n+1)*A_tilde_n(:,:,n) - L_n(:,:,n)'*H*L_n(:,:,n);
        S1(:,n) = q_bold_tilde_n(:,:,n) + A_tilde_n(:,:,n)' * ( S1(:,n+1) + S2(:,:,n+1) * g_tilde_n(:,:,n) ) + ...
            G'*l_n(:,:,n) + L_n(:,:,n)' * (h+H*l_n(:,:,n));
        S0(n) = q_tilde_n(:,:,n) + S0(n+1) + g_tilde_n(:,:,n)'*S1(:,n+1) + (1/2)*g_tilde_n(:,:,n)'*S2(:,:,n+1)*g_tilde_n(:,:,n) ...
            + l_n(:,:,n)' * (h + (1/2)*H*l_n(:,:,n));
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
        u_candidate(:,n) = u_nominal(:,n) + alpha * ( epsilon_n(:,:,n) + Proj_n(:,:,n)*l_n(:,:,n) )...
            + (U_n(:,:,n) + Proj_n(:,:,n)*L_n(:,:,n)) * (x_candidate(:,n) - x_nominal(:,n));  
        k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
        x_mid = x_candidate(:, n) + k1*dt/2;
        k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
        x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
        merit0 = merit0 + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
            + u_candidate(:,n)'*R*u_candidate(:,n) ...
            + sigma*abs( d_n(x_candidate(:, n+1)) );
    end

    %get merit function for alpha = 1
    alpha = 1;
    merit_alpha = 0;
    for n = 1:1:N-1
        disp(n)
        u_candidate(:,n) = u_nominal(:,n) + alpha * ( epsilon_n(:,:,n) + Proj_n(:,:,n)*l_n(:,:,n) )...
            + (U_n(:,:,n) + Proj_n(:,:,n)*L_n(:,:,n)) * (x_candidate(:,n) - x_nominal(:,n));  
        k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
        x_mid = x_candidate(:, n) + k1*dt/2;
        k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
        x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
        merit_alpha = merit_alpha + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
            + u_candidate(:,n)'*R*u_candidate(:,n) ...
            + sigma*abs( d_n(x_candidate(:, n+1)) );
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
                u_candidate(:,n) = u_nominal(:,n) + alpha * ( epsilon_n(:,:,n) + Proj_n(:,:,n)*l_n(:,:,n) )...
                    + (U_n(:,:,n) + Proj_n(:,:,n)*L_n(:,:,n)) * (x_candidate(:,n) - x_nominal(:,n));  
                k1 = x_dot_nonlin(x_candidate(:, n), u_candidate(:, n));
                x_mid = x_candidate(:, n) + k1*dt/2;
                k2 = x_dot_nonlin(x_mid, u_candidate(:, n));
                x_candidate(:, n+1) = x_candidate(:, n) + dt*k2;
                merit_alpha = merit_alpha + x_candidate(:,n)'*Q_n*x_candidate(:,n) ...
                    + u_candidate(:,n)'*R*u_candidate(:,n) ...
                    + sigma*abs( d_n(x_candidate(:, n+1)) );
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