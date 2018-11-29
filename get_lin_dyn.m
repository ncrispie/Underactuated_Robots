function [A, B] = get_lin_dyn(theta)
%get_lin_dyn returns the linearized dynamics: q_dot = A*q + B*u
%   theta: vector of length 4, current position of the four links

L1 = 1; L2 = 1; L3 = 1; L4 = 1; % link lengths
m1 = 1; m2 = 1; m3 = 1; m4 = 1; % link inerias 
b1 = 0.1; b2 = 0.2; b3 = 0.3; b4 = 0.4; % link dampings
g = 9.81; % gravitational constant

%formulate state matrices
A1 = m1*g*L1/2*sin(theta(1)) + m2*g*(L1*sin(theta(1)) + L2/2*sin(theta(1)+theta(2))) ...
    + m3*g*(L1*sin(theta(1))+L2*sin(theta(1)+theta(2))+L3/2*sin(theta(1)+theta(2)+theta(3))) ...
    + m4*g*(L1*sin(theta(1))+L2*sin(theta(1)+theta(2))+L3*sin(theta(1)+theta(2)+theta(3)) + L4/2*sin(theta(1)+theta(2)+theta(3)+theta(4)));

A2 = m2*g*L2/2*sin(theta(1)+theta(2)) + m3*g*(L2*sin(theta(1)+theta(2))+L3*sin(theta(1)+theta(2)+theta(3))) ...
    + m4*(L2*sin(theta(1)+theta(2))+L3*sin(theta(1)+theta(2)+theta(3))+L4/2*sin(theta(1)+theta(2)+theta(3)+theta(4)));

A3 = m3*g*L3/2*sin(theta(1)+theta(2)+theta(3))+m4*g*(L3*sin(theta(1)+theta(2)+theta(3))+L4/2*sin(theta(1)+theta(2)+theta(3)+theta(4)));

A4 = m4*g*L4/2*sin(theta(1)+theta(2)+theta(3)+theta(4));

A_sw_block = [A1/m1 A2/m1 A3/m1 A4/m1; 
              A2/m2 A2/m2 A3/m2 A4/m2;
              A3/m3 A3/m3 A3/m3 A4/m3;
              A4/m4 A4/m4 A4/m4 A4/m4];
          
A_se_block = diag([-b1/m1 -b2/m2 -b3/m3 -b4/m4]);

%state dyanmics in matrix form
A = [zeros(4), eye(4); A_sw_block, A_se_block];

B = [zeros(4,4) ; diag([1/m1 1/m2 1/m3 1/m4])];

C = eye(8);

D = zeros(1,4);

end

