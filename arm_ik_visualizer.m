function theta = arm_ik_visualizer(thetas)
load('dynamics_params.mat');

Theta1=thetas(1);
Theta2=thetas(2);
Theta3 = thetas(3);
Theta4 = thetas(4);


W = eye(4);%weighting matrix (identity all joints penalized equally

%Jacobian (2x3) Partial xee with respect to Theta1,2,3
J= [-L1*sin(Theta1)-L2*sin(Theta1+Theta2)-L3*sin(Theta1+Theta2+Theta3)-L4*sin(Theta1+Theta2+Theta3+Theta4),...
-L2*sin(Theta1+Theta2)-L3*sin(Theta1+Theta2+Theta3)-L4*sin(Theta1+Theta2+Theta3+Theta4),...
-L3*sin(Theta1+Theta2+Theta3)-L4*sin(Theta1+Theta2+Theta3+Theta4),...
-L4*sin(Theta1+Theta2+Theta3+Theta4);...
L1*cos(Theta1)+L2*cos(Theta1+Theta2)+L3*cos(Theta1+Theta2+Theta3)+L4*cos(Theta1+Theta2+Theta3+Theta4),...
L2*cos(Theta1+Theta2)+L3*cos(Theta1+Theta2+Theta3)+L4*cos(Theta1+Theta2+Theta3+Theta4),...
L3*cos(Theta1+Theta2+Theta3)+L4*cos(Theta1+Theta2+Theta3+Theta4),...
L4*cos(Theta1+Theta2+Theta3+Theta4)];
  
%Pseudoinverse given the Jacobian and Weighting matrix
Jpseudo = inv(W)*J'*inv(J*inv(W)*J');

%End - effector location
xo = [L1*cos(Theta1)+L2*cos(Theta1+Theta2)+L3*cos(Theta1+Theta2+Theta3)+L4*cos(Theta1+Theta2+Theta3+Theta4);...
L2*sin(Theta1)+L2*sin(Theta1+Theta2)+L3*sin(Theta1+Theta2+Theta3)+L4*sin(Theta1+Theta2+Theta3+Theta4)];

%Calculate the location of the middle two joints
pointl1 = [L1*cos(Theta1) ; L1*sin(Theta1)];
pointl2 = pointl1 + [L2*cos(Theta1+Theta2);L2*sin(Theta1+Theta2)];
pointl3 = pointl2 + [L3*cos(Theta1+Theta2+Theta3); L3*sin(Theta1+Theta2+Theta3)];

%Plot

axis([-5 5 -5 5])
axis square
line([0,pointl1(1)],[0,pointl1(2,1)],'Color','red')
hold on
line([pointl1(1),pointl2(1)],[pointl1(2,1),pointl2(2,1)],'Color','green')
line([pointl2(1),pointl3(1)],[pointl2(2,1),pointl3(2,1)],'Color','blue')
% line([pointl2(1),xo(1)],[pointl2(2,1),xo(2,1)],'Color','blue')
line([pointl3(1), xo(1)],[pointl3(2,1),xo(2,1)],'Color','magenta')
plot(xo(1),xo(2),'o')
pause(.1)


end