function func_figure1(heter,MC)

means = zeros(1, 22);
mins = zeros(1, 22);
maxs = zeros(1, 22);
Q1s = zeros(1, 22);
Q3s = zeros(1, 22);

% calculate statistic value
for k = 1:22
    data = MC(k,:);
    means(k) = mean(data);        % mean value
    mins(k) = min(data);          % minimum value
    maxs(k) = max(data);          % maximum value
    Q1s(k) = prctile(data, 25);   % first cuartile
    Q3s(k) = prctile(data, 75);   % third cuartile
end

upper_errors = maxs - means;     
lower_errors = means - mins;   

x_vals = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];


figure;
hold on;
set(gcf, 'Color', 'w');

errorbar(x_vals, means, lower_errors, upper_errors, 'o', 'CapSize', 10, 'MarkerSize', 2, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');
plot(x_vals,means,'-ob','LineWidth', 2,'MarkerSize', 2, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');

for k = 1:22

    x_box = [x_vals(k)-0.01, x_vals(k)-0.01, x_vals(k)+0.01, x_vals(k)+0.01];  
    y_box = [Q1s(k), Q3s(k), Q3s(k), Q1s(k)]; 
    
    fill(x_box, y_box, 'c', 'FaceAlpha', 0,  'EdgeColor', 'b', 'LineWidth', 0.6);

end

percentages = [5, 10, 15, 20];  
title(['Battery Capacity Variation ', num2str(percentages(heter)), '%']);

xlabel('Normalized Aggregate Converter Rating');
ylabel('Battery Power Utilization(%)');

xticks(x_vals);
xticklabels({'0.02','0.04','0.06','0.08','0.1','0.12','0.14','0.16','0.18','0.2','0.22', '0.24','0.26', '0.28','0.3', '0.4','0.5', '0.6','0.7', '0.8', '0.9', '1.0'});
xtickangle(90);

xlim([0 1.02])
ylim([35 100])

grid on;

hold off;