%pick cost constants
Q = eye(8);
R = eye(4);

% initial state
theta = [pi/4 ; pi/6 ; -pi/3 ; pi/4];

%TODO" Compute cost

[A, B] = get_lin_dyn(theta);

%this yields
qdot = A*q + B*u;
%q =
%[theta1,theta2,theta3,theta3,theta4,...
%theta1_dot,theta2_dot,theta3_dot,theta4_dot]'
%u = [T1 T2 T3 T4]'

%TODO
%write out constraint matrix


