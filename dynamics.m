load('dynamics_params.mat');

% pick cost constants
Q = eye(8);
R = eye(4);

% initial state
% theta = [pi/4 ; pi/6 ; -pi/3 ; pi/4];
% theta = 0 is straight down. thetapositive is counter clockwise.
x_0 = [pi/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
u_0 = [0 ; 0 ; 0 ; 0];

%TODO" Compute cost

% [A, B] = get_lin_dyn(x);

% this yields
% x_dot_lin = A*x + B*u

% x_dot_nonlin = get_dyn(x, u)

% get symbolic derivatives and make lambda function accessors
[x, u, x_dot_nonlin_sym] = get_dyn();
x_dot_nonlin = @(x_eval, u_eval) double(subs(x_dot_nonlin_sym,  [x ; u], [x_eval ; u_eval]));
x_dot_lin = @(x_about, u_about, x_eval, u_eval) ...
    double(subs(taylor(x_dot_nonlin_sym, [x ; u], [x_about ; u_about],'Order', 2), ...
    [x ; u], [x_eval ; u_eval]));

% usage
x_dot_lin(x_0, u_0, x_0, u_0);
x_dot_nonlin_0 = x_dot_nonlin(x_0, u_0);


% linearized constraints:
% remember, linearize about the nominal trajectory
% apply to C_N and d_N also
C_n = @(x) [L1*sin(x(1))+L2*sin(x(1)+x(2))+L3*sin(x(1)+x(2)+x(3))+L4*sin(x(1)+x(2)+x(3)+x(4)) ...
            L2*sin(x(1)+x(2))+L3*sin(x(1)+x(2)+x(3))+L4*sin(x(1)+x(2)+x(3)+x(4)) ...
            L3*sin(x(1)+x(2)+x(3))+L4*sin(x(1)+x(2)+x(3)+x(4)) ...
            L4*sin(x(1)+x(2)+x(3)+x(4)) ; ...
            -L1*cos(x(1))-L2*cos(x(1)+x(2))-L3*cos(x(1)+x(2)+x(3))-L4*cos(x(1)+x(2)+x(3)+x(4)) ...
            -L2*cos(x(1)+x(2))-L3*cos(x(1)+x(2)+x(3))-L4*cos(x(1)+x(2)+x(3)+x(4)) ...
            -L3*cos(x(1)+x(2)+x(3))-L4*cos(x(1)+x(2)+x(3)+x(4)) ...
            -L4*cos(x(1)+x(2)+x(3)+x(4))];
d_n = [0 ; 0];



