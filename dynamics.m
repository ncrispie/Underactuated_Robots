load('dynamics_params.mat');

N = 51; % number of points, meaning N-1 windows

% pick cost constants
% we're leaving our various other costs: 
% q_N, q_bold_N, q_n, q_bold_n, r_n, P_n
Q_n = diag([1 1 1 1 1 1 1 1]);
Q_N = Q_n;
R = diag([1 1 1 1]);

% initial state
% theta = 0 is straight down. theta positive is counter clockwise.
x_0 = [pi/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
u_0 = [0 ; 0 ; 0 ; 0];
x_1 = [pi/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
u_1 = [0 ; 0 ; 0 ; 0];

% get symbolic derivatives
[x, u, x_dot_nonlin_sym] = get_dyn();
A_sym = jacobian(x_dot_nonlin_sym, x);
B_sym = jacobian(x_dot_nonlin_sym, u);

% lambda accessor functions
x_dot_nonlin = @(x_eval, u_eval) double(subs(x_dot_nonlin_sym,  [x ; u], [x_eval ; u_eval]));
x_dot_lin = @(x_about, u_about, x_eval, u_eval) ...
    double(subs(taylor(x_dot_nonlin_sym, [x ; u], [x_about ; u_about],'Order', 2), ...
    [x ; u], [x_eval ; u_eval]));
A = @(x_about, u_about) double(subs(A_sym, [x ; u], [x_about ; u_about]));
B = @(x_about, u_about) double(subs(B_sym, [x ; u], [x_about ; u_about]));

% usage
% remember, linearize about the nominal trajectory
x_dot_lin(x_0, u_0, x_0, u_0);
x_dot_nonlin_0 = x_dot_nonlin(x_0, u_0);
A_0 = A(x_0, u_0);
B_0 = B(x_0, u_0);

% nonlinear constraint
% g2(x_n,n) = 0
% this also applies at N! g3(x_N,N) = g2(x_n,n)
% ee is end effector
ee_start = [1 ; -1];
ee_end = [2 ; 1];
g2 = @(x_n, n) [ee_start(1) + (ee_end(1) - ee_start(1))*(n-1)/(N-1) - ...
    (L1*cos(x_n(1)) + L2*cos(x_n(1)+x_n(2)) + L3*cos(x_n(1)+x_n(2)+x_n(3)) ...
    + L4*cos(x_n(1)+x_n(2)+x_n(3)+x_n(4))) ; ...
    ee_start(2) + (ee_end(2) - ee_start(2))*(n-1)/(N-1) - ...
    L1*sin(x_n(1)) + L2*sin(x_n(1)+x_n(2)) + L3*sin(x_n(1)+x_n(2)+x_n(3)) ...
    + L4*sin(x_n(1)+x_n(2)+x_n(3)+x_n(4))];

% linearized constraint
% remember to linearize about the nominal trajectory!
% the x_n passed in here should come from the nominal trajectory
% apply to C_N and d_N also
C_n = @(x_n) [L1*sin(x_n(1))+L2*sin(x_n(1)+x_n(2))+L3*sin(x_n(1)+x_n(2)+x_n(3))+L4*sin(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
              L2*sin(x_n(1)+x_n(2))+L3*sin(x_n(1)+x_n(2)+x_n(3))+L4*sin(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
              L3*sin(x_n(1)+x_n(2)+x_n(3))+L4*sin(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
              L4*sin(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
              0 0 0 0 ; ...
             -L1*cos(x_n(1))-L2*cos(x_n(1)+x_n(2))-L3*cos(x_n(1)+x_n(2)+x_n(3))-L4*cos(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
             -L2*cos(x_n(1)+x_n(2))-L3*cos(x_n(1)+x_n(2)+x_n(3))-L4*cos(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
             -L3*cos(x_n(1)+x_n(2)+x_n(3))-L4*cos(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
             -L4*cos(x_n(1)+x_n(2)+x_n(3)+x_n(4)) ...
              0 0 0 0];
d_n = [0 ; 0];

M_n = C_n(x_1) * A(x_1,u_1) * B(x_0,u_0);
N_n = C_n(x_1) * A(x_1,u_1) * A(x_0,u_0);

