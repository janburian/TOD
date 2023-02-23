function [] = plot_ii(trajectory, measurement, kalman, kalman_const, deterministic_estimator, kalman_inovace, kalman_const_inovace)

%% POROVN�N� ODHADU POLOHY
figure
hold on
plot(0:length(kalman)-1, trajectory(1, :), 'LineWidth', 1.3)
plot(0:length(kalman)-1, [NaN measurement(2:end)], 'o', 'LineWidth', 1.3)
plot(0:length(kalman)-1, kalman(1, :))
plot(0:length(kalman)-1, kalman_const(1, :))
plot(0:length(kalman)-1, deterministic_estimator(1, :))
title('Porovn�n� odhadu polohy')
legend('Skute�n� trajektorie', 'M��en�', 'Kalman�v filtr', 'Ust�len� Kalman�v filtr', 'Deterministick� rekonstruktor')
xlabel('�as')
ylabel('Poloha')

%% CHYBY ODHAD� POLOHY
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
title('Porovn�n� chyby odhadu polohy')
legend([p1, p2, p3], 'Kalman�v filtr', 'Ust�len� Kalman�v filtr', 'Deterministick� rekonstruktor')
xlabel('�as')
ylabel('Chyba odhadu polohy')

%% POROVN�N� ODHADU RYCHLOSTI
figure
hold on
plot(0:length(kalman)-1, trajectory(2, :), 'LineWidth', 1.3)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
plot(0:length(kalman)-1, kalman(2, :))
plot(0:length(kalman)-1, kalman_const(2, :))
plot(0:length(kalman)-1, deterministic_estimator(2, :))
title('Porovn�n� odhadu rychlosti')
legend('Skute�n� trajektorie', 'Kalman�v filtr', 'Ust�len� Kalman�v filtr', 'Deterministick� rekonstruktor')
xlabel('�as')
ylabel('Rychlost')

%% CHYBY ODHAD� RYCHLOSTI
figure
hold on
plot(0:length(kalman)-1, zeros(1, length(measurement)), ':');
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
p1 = plot(0:length(kalman)-1, e_kalman(2, :));
p2 = plot(0:length(kalman)-1, e_kalman_const(2, :));
p3 = plot(0:length(kalman)-1, e_deterministic_estimator(2, :));
title('Porovn�n� chyby odhadu rychlosti')
legend([p1, p2, p3], 'Kalman�v filtr', 'Ust�len� Kalman�v filtr', 'Deterministick� rekonstruktor')
xlabel('�as')
ylabel('Chyba odhadu polohy')

%% INOVA�N� POSLOUPNOSTI
figure
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
plot(0:length(kalman)-1, kalman_inovace);
plot(0:length(kalman)-1, kalman_const_inovace);
title('Porovn�n� inova�n�ch posloupnost�')
legend('Kalman�v filtr', 'Ust�len� Kalman�v filtr')
xlabel('�as')
ylabel('Inovace')

end

