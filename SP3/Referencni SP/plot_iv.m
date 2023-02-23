function [] = plot_iv(mean_e_kalman_t,mean_e_kalman_sim, mean_e_kalman_const_t, mean_e_kalman_const_sim, mean_e_det_est_t, mean_e_det_est_sim)

%% POLOHA
figure
hold on
% ax = gca;
% ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
% %mean_e_kalman_t = zeros(2,;
% plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_sim(1, :), ':', 'LineWidth', 1.5)
% ax = gca;
% ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
% plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_t(1, 1:2:length(mean_e_kalman_t)))

mean_e_kalman_const_t = cell2mat(mean_e_kalman_const_t);
plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_const_sim(1, :), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_const_t(1, 1:2:length(mean_e_kalman_const_t)))

mean_e_det_est_t = cell2mat(mean_e_det_est_t);
plot(0:length(mean_e_kalman_sim)-1, mean_e_det_est_sim(1, :), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, mean_e_det_est_t(1, 1:2:length(mean_e_det_est_t)))
title('Èasový prùbìh støední kvadratické chyby (Poloha)')
xlabel('Èas')
ylabel('Støední kvadratická chyba polohy')
legend('Kalman - simulace', 'Kalman - teorie', 'Ustálený Kalman - simulace', 'Ustálený Kalman - teorie', 'Deterministický rekonstruktor - simulace', 'Deterministický rekonstruktor - teorie')

%% RYCHLOST
figure
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_sim(2, :), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_t(2, 2:2:length(mean_e_kalman_t)))

plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_const_sim(2, :), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, mean_e_kalman_const_t(2, 2:2:length(mean_e_kalman_const_t)))

plot(0:length(mean_e_kalman_sim)-1, mean_e_det_est_sim(2, :), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, mean_e_det_est_t(2, 2:2:length(mean_e_det_est_t)))
title('Èasový prùbìh støední kvadratické chyby (Rychlost)')
xlabel('Èas')
ylabel('Støední kvadratická chyba rychlosti')
legend('Kalman - simulace', 'Kalman - teorie', 'Ustálený Kalman - simulace', 'Ustálený Kalman - teorie', 'Deterministický rekonstruktor - simulace', 'Deterministický rekonstruktor - teorie')


%% STOPA
a_mean_e_kalman_t = mean_e_kalman_t(1, :);
b_mean_e_kalman_t = mean_e_kalman_t(2, :);
mean_e_kalman_t = [a_mean_e_kalman_t(1:2:length(a_mean_e_kalman_t));
                   b_mean_e_kalman_t(2:2:length(b_mean_e_kalman_t))];
a_mean_e_kalman_const_t = mean_e_kalman_const_t(1, :);
b_mean_e_kalman_const_t = mean_e_kalman_const_t(2, :);
mean_e_kalman_const_t = [a_mean_e_kalman_const_t(1:2:length(a_mean_e_kalman_const_t));
                   b_mean_e_kalman_const_t(2:2:length(b_mean_e_kalman_const_t))];
a_mean_e_det_est_t = mean_e_det_est_t(1, :);
b_mean_e_det_est_t = mean_e_det_est_t(2, :);
mean_e_det_est_t = [a_mean_e_det_est_t(1:2:length(a_mean_e_det_est_t));
                   b_mean_e_det_est_t(2:2:length(b_mean_e_det_est_t))];               


figure
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex + 2;
plot(0:length(mean_e_kalman_sim)-1, sum(mean_e_kalman_sim), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, sum(mean_e_kalman_t))

plot(0:length(mean_e_kalman_sim)-1,sum(mean_e_kalman_const_sim), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, sum(mean_e_kalman_const_t))

plot(0:length(mean_e_kalman_sim)-1, sum(mean_e_det_est_sim), ':', 'LineWidth', 1.5)
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(0:length(mean_e_kalman_sim)-1, sum(mean_e_det_est_t))
title('Èasový prùbìh støední kvadratické chyby')
xlabel('Èas')
ylabel('Støední kvadratická chyba')
legend('Kalman - simulace', 'Kalman - teorie', 'Ustálený Kalman - simulace', 'Ustálený Kalman - teorie', 'Deterministický rekonstruktor - simulace', 'Deterministický rekonstruktor - teorie')


end

