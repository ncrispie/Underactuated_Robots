function [x, u, x_dot] = get_dyn()
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
A_coeffs = sym('A_coeffs', [3 3]);
A_coeffs(1,1) = (I1 + (m1/4 + m2 + m3)*L1^2);
A_coeffs(1,2) = ((m2/2 + m3)*L1*L2);
A_coeffs(1,3) = ((m3/2)*L1*L3);
A_coeffs(2,1) = (A_coeffs(1,2));
A_coeffs(2,2) = (I2 + (m2/4 + m3)*L2^2);
A_coeffs(2,3) = ((m3/2)*L2*L3);
A_coeffs(3,1) = (A_coeffs(1,3));
A_coeffs(3,2) = (A_coeffs(2,3));
A_coeffs(3,3) = (I3 + (m3/4)*L3^2);
A_state = sym('A_state', [3 3]);
for i = 1:3
    for j = 1:3
        A_state(i,j) = cos(theta(i) - theta(j));
    end
end
A = A_state .* A_coeffs;

% B_coeffs = zeros(4,4);
B_coeffs = sym('B_coeffs', [3 3]);
B_coeffs(1,1) = (-b1-b2);
B_coeffs(1,2) = (b2+(m2/2 + m3)*L1*L2);
B_coeffs(1,3) = ((m3/2)*L1*L3);
B_coeffs(2,1) = (b2-(m2/2 + m3)*L1*L2);
B_coeffs(2,2) = (-b2-b3);
B_coeffs(2,3) = (b3+(m3/2)*L2*L3);
B_coeffs(3,1) = (-B_coeffs(1,3));
B_coeffs(3,2) = (b3-(m3/2)*L2*L3);
B_coeffs(3,3) = (-b3);
B_state = sym('B_state', [3 3]);
for i = 1:3
    for j = 1:3
        if i ~= j
            k = max(i,j);
            n = min(i,j);
            B_state(i,j) = theta_dot(j) * sin(theta(k) - theta(n));
        else
            B_state(i,j) = 1;
        end
    end
end
B = B_state .* B_coeffs;

% C_coeffs = zeros(4,1);
C_coeffs = sym('C_coeffs', [3 1]);
C_coeffs(1) = g*L1*((1/2)*m1+m2+m3);
C_coeffs(2) = g*L2*((1/2)*m2+m3);
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
theta_ddot = inv(A) * -(B*theta_dot + C + D);

x_dot = [theta_dot ; theta_ddot];

end

