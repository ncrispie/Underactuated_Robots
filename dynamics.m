%initial conditions
L1 = 1; L2 = 1; L3 = 1; L4 = 1;
m1 = 1; m2 = 1; m3 = 1; m4 = 1;
g = 9.81;

%pick cost constants
Q = eye(8);
R = eye(4);
%TODO" Compute cost

%formulate state matrices
A1 = m1*g*L1/2*sin(theta1) + m2*g*(L1*sin(theta1) + L2/2*sin(theta1+theta2)) ...
    + m3*g*(L1*sin(theta1)+L2*sin(theta1+theta2)+L3/2*sin(theta1+theta2+theta3)) ...
    + m4*g*(L1*sin(theta1)+L2*sin(theta1+theta2)+L3*sin(theta1+theta2+theta3) + L4/2*sin(theta1+theta2+theta3+theta4));

A2 = m2*g*L2/2*sin(theta1+theta2) + m3*g*(L2*sin(theta1+theta2)+L3*sin(theta1+theta2+theta3)) ...
    + m4*(L2*sin(theta1+theta2)+L3*sin(theta1+theta2+theta3)+L4/2*sin(theta1+theta2+theta3+theta4));

A3 = m3*g*L3/2*sin(theta1+theta2+theta3)+m4*g*(L3*sin(theta1+theta2+theta3)+L4/2*sin(theta1+theta2+theta3+theta4));

A4 = m4*g*L4/2*sin(theta1+theta2+theta3+theta4);

A_sw_block = [A1/m1 A2/m1 A3/m1 A4/m1; 
              A2/m2 A2/m2 A3/m2 A4/m2;
              A3/m3 A3/m3 A3/m3 A4/m3;
              A4/m4 A4/m4 A4/m4 A4/m4];
          
A_se_block = [-b/m1 -b/m2 -b/m3 -b/m4].*eye(4);

%state dyanmics in matrix form
A = [zeros(4), eye(4); A_sw_block, A_se_block];

B = [0 0 0 0 1/m1 1/m2 1/m3 1/m4]';

C = eye(8);

D = zeros(1,4);

%this yields
qdot = A*q + B*u;
%q =
%[theta1,theta2,theta3,theta3,theta4,...
%theta1_dot,theta2_dot,theta3_dot,theta4_dot]'
%u = [T1 T2 T3 T4]'

%TODO
%write out constraint matrix



