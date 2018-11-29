function [theta_ddot] = get_dyn(theta, theta_dot, u)
%get_dyn returns the theta_dot based on nonlinear dynamics
%   theta: vector of length 4, current position of the four links
%   u: vector of length 4, input at each of the four links

% construct A*x_ddot = B*x_dot + C*X + D

L1 = 1; L2 = 1; L3 = 1; L4 = 1; % link lengths
m1 = 1; m2 = 1; m3 = 1; m4 = 1; % link masses 
I1 = m1 * L1^2 / 12; I2 = m2 * L2^2 / 12; I3 = m3 * L3^2 / 12; I4 = m4 * L4^2 / 12; % link inertias
b1 = 0.1; b2 = 0.2; b3 = 0.3; b4 = 0.4; % link dampings
g = 9.81; % gravitational constant

A = zeros(4,4);
A(1,1) = (I1 + (m1/4 + m2 + m3 + m4)*L1^2);
A(1,2) = ((m2/2 + m3 + m4)*L1*L2);
A(1,3) = ((m3/2 + m4)*L1*L3);
A(1,4) = (m4*L1*L4/2);
A(2,1) = (A(1,2));
A(2,2) = (I2 + (m2/4 + m3 + m4)*L2^2);
A(2,3) = ((m3/2 + m4)*L2*L3);
A(2,4) = (m4*L2*L4/2);
A(3,1) = (A(1,3));
A(3,2) = (A(2,3));
A(3,3) = (I3 + (m3/4 + m4)*L3^2);
A(3,4) = (m4*L3*L4/2);
A(4,1) = (A(1,4));
A(4,2) = (A(2,4));
A(4,3) = (A(3,4));
A(4,4) = (I4 + m4*L4^2/2);
for i = 1:4
    for j = 1:4
        A(i,j) = A(i,j) * cos(theta(i) - theta(j));
    end
end

B = zeros(4,4);
B(1,1) = (-b1-b2);
B(1,2) = (b2+(m2/2 + m3+m4)*L1*L2);
B(1,3) = ((m3/2 + m4)*L1*L2);
B(1,4) = (m4*L1*L4/2);
B(2,1) = (b2-(m2/2 + m3 + m4)*L1*L2);
B(2,2) = (-b2-b3);
B(2,3) = (b3+(m3/2 + m4)*L2*L3);
B(2,4) = (m4*L2*L4/2);
B(3,1) = (-B(1,3));
B(3,2) = (b3-(m3/2 + m4)*L2*L3);
B(3,3) = (-b3-b4);
B(3,4) = (b4+m4*L3*L4/2);
B(4,1) = (-B(1,4));
B(4,2) = (-B(2,4));
B(4,3) = (b4-m4*L3*L4/2);
B(4,4) = (-b4);
for i = 1:4
    for j = 1:4
        if i ~= j
            k = max(i,j);
            n = min(i,j);
            A(i,j) = A(i,j) * theta_dot(j) * sin(theta(k) - theta(n));
        end
    end
end

C = zeros(4,1);
C(1) = g*L1*((1/2)*m1+m2+m3+m4);
C(2) = g*L2*((1/2)*m2+m3+m4);
C(3) = g*L3*((1/2)*m3+m4);
C(4) = g*L4*((1/2)*m4);
for i = 1:4
    C(i) = C(i) * sin(theta(i));
end

D = u;

theta_ddot = inv(A) * (B*theta_dot + C + D);

end

