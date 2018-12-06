function plot_joints(thetas)
%   thetas should be 3xN
figure(1)
hold on
dt = 0.1;
plot(dt*(0:size(thetas,2)-1), thetas(1,:)', 'b-', 'LineWidth', 2);
plot(dt*(0:size(thetas,2)-1), thetas(2,:)', 'r--', 'LineWidth', 2);
plot(dt*(0:size(thetas,2)-1), thetas(3,:)', 'g--', 'LineWidth', 2);
xlabel('Time')
ylabel('Joint Angles')
legend({'Joint 1', 'Joint 2', 'Joint 3'}, 'Location', 'best')
end