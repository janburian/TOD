function [] = plot_ii(trajectory, measurement, kalman, kalman_const, deterministic_estimator, kalman_inovace, kalman_const_inovace)

%% POROVNÁNÍ ODHADU POLOHY
figure
hold on
plot(0:length(kalman)-1, trajectory(1, :), 'LineWidth', 1.3)
plot(0:length(kalman)-1, [NaN measurement(2:end)], 'o', 'LineWidth', 1.3)
plot(0:length(kalman)-1, kalman(1, :))
plot(0:length(kalman)-1, kalman_const(1, :))
plot(0:length(kalman)-1, deterministic_estimator(1, :))
title('Porovnání odhadu polohy')
legend('Skuteèná trajektorie', 'Mìøení', 'Kalmanùv filtr', 'Ustálený Kalmanùv filtr', 'Deterministický rekonstruktor')
xlabel('Èas')
ylabel('Poloha')

%% CHYBY ODHADÙ POLOHY
e_kalman = trajectory - kalman;
e_kalman_const = trajectory - kalman_const;
e_deterministic_estimator = trajectory - deterministic_estimator;

figure
hold on
plot(0:length(kalman)-1, zeros(1, length(measurement)), ':');
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
p1 = plot(0:length(kalman)-1, e_kalman(1, :));
p2 = plot(0:length(kalman)-1, e_kalman_const(1, :));
p3 = plot(0:length(kalman)-1, e_deterministic_estimator(1, :));
title('Porovnání chyby odhadu polohy')
legend([p1, p2, p3], 'Kalmanùv filtr', 'Ustálený Kalmanùv filtr', 'Deterministický rekonstruktor')
xlabel('Èas')
ylabel('Chyba odhadu polohy')

%% POROVNÁNÍ ODHADU RYCHLOSTI
figure
hold on
plot(0:length(kalman)-1, trajectory(2, :), 'LineWidth', 1.3)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
plot(0:length(kalman)-1, kalman(2, :))
plot(0:length(kalman)-1, kalman_const(2, :))
plot(0:length(kalman)-1, deterministic_estimator(2, :))
title('Porovnání odhadu rychlosti')
legend('Skuteèná trajektorie', 'Kalmanùv filtr', 'Ustálený Kalmanùv filtr', 'Deterministický rekonstruktor')
xlabel('Èas')
ylabel('Rychlost')

%% CHYBY ODHADÙ RYCHLOSTI
figure
hold on
plot(0:length(kalman)-1, zeros(1, length(measurement)), ':');
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
p1 = plot(0:length(kalman)-1, e_kalman(2, :));
p2 = plot(0:length(kalman)-1, e_kalman_const(2, :));
p3 = plot(0:length(kalman)-1, e_deterministic_estimator(2, :));
title('Porovnání chyby odhadu rychlosti')
legend([p1, p2, p3], 'Kalmanùv filtr', 'Ustálený Kalmanùv filtr', 'Deterministický rekonstruktor')
xlabel('Èas')
ylabel('Chyba odhadu polohy')

%% INOVAÈNÍ POSLOUPNOSTI
figure
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
plot(0:length(kalman)-1, kalman_inovace);
plot(0:length(kalman)-1, kalman_const_inovace);
title('Porovnání inovaèních posloupností')
legend('Kalmanùv filtr', 'Ustálený Kalmanùv filtr')
xlabel('Èas')
ylabel('Inovace')

end

