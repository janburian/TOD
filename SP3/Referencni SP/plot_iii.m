function [] = plot_iii(esimate_var_kalman_inovace, esimate_var_kalman_const_inovace, x_axis)

len = length(esimate_var_kalman_inovace);
if len ~= length(esimate_var_kalman_const_inovace)
    disp('ERR: Vectors must have same length.')
end

x_axis = x_axis - ones(1, length(x_axis));

figure
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
plot(x_axis, esimate_var_kalman_inovace, 'o')
plot(x_axis, esimate_var_kalman_const_inovace, 'x')
title('Porovnání odhadu variancí inovací')
legend('Kalmanùv filtr', 'Ustálený Kalmanùv filtr')
xlabel('Èas')
ylabel('Variance inovace')

end

