function [x, u, x_dot] = get_dyn2()
%get_dyns returns the x_dot based on nonlinear dynamics
%   x: vector of length 8, current positions & velocities of the four links
%   x = [t1 ; t2 ; t3 ; t4 ; tdot1 ; tdot2 ; tdot3 ; tdot4]
%   u: vector of length 4, input at each of the four links

% construct A*x_ddot = B*x_dot + C + D
% refer to Dynamic Model and Motion Control Analysis of Three-link 
% Gymnastic Robot on Horizontal Bar
% we basically pattern matching the coefficients from this paper

syms x1 x2 x3 x4 x5 x6
x = [x1 ; x2 ; x3 ; x4 ; x5 ; x6];
theta = [x1 ; x2 ; x3];
theta_dot = [x4 ; x5 ; x6];

syms u1 u2 u3
u = [u1 ; u2 ; u3];

syms b1 b2 b3
syms g
syms I1 I2 I3
syms L1 L2 L3
syms m1 m2 m3

% load('dynamics_params.mat');

% A_coeffs = zeros(4,4);
A = diag([I1 I2 I3]);

% B_coeffs = zeros(4,4);
B = diag([-b1 -b2 -b3]);

% C_coeffs = zeros(4,1);
C_coeffs = sym('C_coeffs', [3 1]);
C_coeffs(1) = g*L1*((1/2)*m1);
C_coeffs(2) = g*L2*((1/2)*m2);
C_coeffs(3) = g*L3*((1/2)*m3);
C_state = sym('C_state', [3 1]);
for i = 1:3
    % note that negative sin is a difference from the paper
    % in the paper, zero theta is up and positive theta is clockwise
    % we're adding pi to our thetas so that theta is now down, but positive
    % is still clockwise
    % this is the only effect
    C_state(i) = sin(theta(i));
end
C = C_state .* C_coeffs;

D = u;

% sub in params from file
load('dynamics_params.mat');
A = subs(A); B = subs(B); C = subs(C);

% A\b is the same as inv(A)*b
% note this is different from A/b
theta_ddot = inv(A) * (B*theta_dot + C + D);

x_dot = [theta_dot ; theta_ddot];

end

