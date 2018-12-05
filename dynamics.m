load('dynamics_params.mat');

%pick cost constants
Q = eye(8);
R = eye(4);

% initial state
% theta = [pi/4 ; pi/6 ; -pi/3 ; pi/4];
% theta = 0 is straight down. thetapositive is counter clockwise.
q_0 = [pi/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
u_0 = [0 ; 0 ; 0 ; 0];

%TODO" Compute cost

% [A, B] = get_lin_dyn(q);

% this yields
% q_dot_lin = A*q + B*u

% q_dot_nonlin = get_dyn(q, u)

% get symbolic derivatives and make lambda function accessors
[q, u, q_dot_nonlin_sym] = get_dyn();
q_dot_nonlin = @(q_eval, u_eval) double(subs(q_dot_nonlin_sym,  [q ; u], [q_eval ; u_eval]));
q_dot_lin = @(q_about, u_about, q_eval, u_eval) ...
    double(subs(taylor(q_dot_nonlin_sym, [q ; u], [q_about ; u_about],'Order', 2), ...
    [q ; u], [q_eval ; u_eval]));

% usage
q_dot_lin(q_0, u_0, q_0, u_0);
q_dot_nonlin_0 = q_dot_nonlin(q_0, u_0);


% linearized constraints:
% remember, linearize about the nominal trajectory
% apply to C_N and d_N also
C_n = @(q) [L1*sin(q(1))+L2*sin(q(1)+q(2))+L3*sin(q(1)+q(2)+q(3))+L4*sin(q(1)+q(2)+q(3)+q(4)) ...
            L2*sin(q(1)+q(2))+L3*sin(q(1)+q(2)+q(3))+L4*sin(q(1)+q(2)+q(3)+q(4)) ...
            L3*sin(q(1)+q(2)+q(3))+L4*sin(q(1)+q(2)+q(3)+q(4)) ...
            L4*sin(q(1)+q(2)+q(3)+q(4)) ; ...
            -L1*cos(q(1))-L2*cos(q(1)+q(2))-L3*cos(q(1)+q(2)+q(3))-L4*cos(q(1)+q(2)+q(3)+q(4)) ...
            -L2*cos(q(1)+q(2))-L3*cos(q(1)+q(2)+q(3))-L4*cos(q(1)+q(2)+q(3)+q(4)) ...
            -L3*cos(q(1)+q(2)+q(3))-L4*cos(q(1)+q(2)+q(3)+q(4)) ...
            -L4*cos(q(1)+q(2)+q(3)+q(4))];
d_n = [0 ; 0];



