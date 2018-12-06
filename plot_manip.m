function plot_manip(theta)
figure(1)
hold on
xlim([-3 3])
ylim([-3 3])
xlabel('x')
ylabel('y')
view([90,-90])
cur_pt = [0 ; 0];
for i = 1:size(theta,1)
    next_pt = cur_pt + [cos(theta(i)) ; sin(theta(i))];
    p = plot(cur_pt(1), cur_pt(2), 'k.', 'MarkerSize', 18);
    p.HandleVisibility = 'off';
    l = line([cur_pt(1) next_pt(1)], [cur_pt(2) next_pt(2)],'Color','red' ...
                                                       , 'LineWidth',2);
    if i ~= 1
        l.HandleVisibility = 'off';
    end
    cur_pt = next_pt;
end
end

% plot_manip(x_nominal_history(1:3,end,1))
% plot_manip(x_nominal_history(1:3,end,5))
% plot_manip(x_nominal_history(1:3,end,10))
% legend({'Initial Trajectory', 'Iteration 5', 'Iteration 10 (final)'}, 'Location', 'northeast')