load('final_results.mat');
plot(0:10, reshape(sum(sum(abs(u_nominal_history))), [11 1]), 'k.', 'MarkerSize', 18);
xlabel('Iteration');
ylabel('Total Absolute Control');