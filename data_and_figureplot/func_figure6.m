function func_figure6(x,y)

x1 = round(x./2700, 2); 
y1 = round(y./2700, 2); 

figure;
set(gcf, 'Color', 'w');

plot(x1, y1, 'k-', 'LineWidth', 2);
axis([-3 3 -2 2]);  
grid on;
hold on;

title('Normalized Aggregate Converter Rating = 0.5');
xlabel('Reactive Power (p.u.)');
ylabel('Active Power (p.u.)');

axis equal;

% plot(0, 0, 'ko', 'MarkerSize', 1, 'MarkerFaceColor', 'k');

for i = 1:length(x1)
    plot([x1(i), x1(i)], [y1(i), 0], 'k--'); 
    plot([x1(i), 0], [y1(i), y1(i)], 'k--');  
end

for i = 1:length(x1)
    if x1(i) ~= 0
    text(x1(i), -0.2, num2str(x1(i)), 'HorizontalAlignment', 'center');
    end
end

for i = 1:length(y1)
    if y1(i) ~= 0
    text(-0.1, y1(i)-0.2, num2str(y1(i)), 'HorizontalAlignment', 'right');
    end
end

hold off;
end

