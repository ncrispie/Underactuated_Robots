load('dynamics_params.mat');

%pick cost constants
Q = eye(8);
R = eye(4);

% initial state
% theta = [pi/4 ; pi/6 ; -pi/3 ; pi/4];
% theta = 0 is straight down. thetapositive is counter clockwise.
q = [pi/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
u = [0 ; 0 ; 0 ; 0];

%TODO" Compute cost

% [A, B] = get_lin_dyn(q);

% this yields
% q_dot_lin = A*q + B*u

q_dot_nonlin = get_dyn(q, u)

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



