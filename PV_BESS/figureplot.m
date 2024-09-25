x = linspace(1,33,33);

figure;
hold on;
set(gcf, 'Color', 'w');

% for i=1:33
%     datapro(2,i) = sqrt(datapro(1,i));
% end

num_colors = 4;
linecolors = [linspace(0.9, 0.1, num_colors)', linspace(0.8, 0.2, num_colors)',  linspace(1,0.7, num_colors)'];

plot(x, datapro(1,:), '-o', 'MarkerFaceColor', [0.4	0.8	1], 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '1 WT + 1 2BESS');
plot(x, datapro(2,:), '-o', 'MarkerFaceColor', linecolors(2, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '2 WT + 2 2BESS');
plot(x, datapro(3,:), '-o', 'MarkerFaceColor', linecolors(3, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '3 WT + 3 2BESS');
plot(x, datapro(4,:), '-o', 'MarkerFaceColor', linecolors(4, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '4 WT + 4 2BESS');
% plot(x, datapro(5,:), '-o', 'MarkerFaceColor', linecolors(5, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '4 BESS');
% plot(x, datapro(6,:), '-o', 'MarkerFaceColor', linecolors(6, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '5 BESS');
% plot(x, datapro(7,:), '-o', 'MarkerFaceColor', linecolors(7, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '6 BESS');
% plot(x, datapro(8,:), '-o', 'MarkerFaceColor', linecolors(8, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '7 BESS');
% plot(x, datapro(9,:), '-o', 'MarkerFaceColor', linecolors(9, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '8 BESS');
% plot(x, datapro(10,:), '-o', 'MarkerFaceColor', linecolors(10, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', '9 BESS');


legend('show');

xlim([1 35]);
ylim([0.92 1]); 

title('Buses voltage profile');
xlabel('Bus Number');
ylabel('Voltage');
% ylabel('System Power Efficiency(%)');

grid on;

hold off;