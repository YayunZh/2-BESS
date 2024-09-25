function func_figure2(heter,mac)

load("dcacUp.mat");

x = [0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];

figure;
hold on;
set(gcf, 'Color', 'w');

plot(x, dcacUp(heter, :), '-om', 'LineWidth', 2, 'DisplayName', 'DC-AC');
plot(x, mac, '-ob', 'LineWidth', 2, 'DisplayName', 'MAC-2BESS');

legend('show');

xlim([0 1.05]);
ylim([80 100]); 

percentages = [5, 10, 15, 20];  
title(['Battery Capacity Variation ', num2str(percentages(heter)), '%']);

% title('Battery Capacity Variation 20%');
xlabel('Normalized Aggregate Converter Rating');
ylabel('Battery Power Utilization(%)');

grid on;

hold off;

end

