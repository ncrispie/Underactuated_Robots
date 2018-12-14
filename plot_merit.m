% plot merit
% load('final_results.mat');
plot(0:10, merit_history, 'k.', 'MarkerSize', 18);
xlabel('Iteration');
ylabel('Total Trajectory Cost');