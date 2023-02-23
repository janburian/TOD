%% 3. SEMESTRÁLNÍ PRÁCE Z PØEDMÌTU KKY/TOD
% autor: Tomáš Honzík
% datum: 15.12.2019

%% ZADÁNÍ
close all; clear; clc;

T = 1; % Èasová konstanta

A = [1 T; % Matice dynamiky
     0 1];

C = [1 0]; % Matice výstupu

G = eye(2); % Matice šumu dynamiky

% Momenty poèáteèní podmínky
E_x0 = [10 1]';
P_x0 = [10 0;
        0  1];

% Momenty šumu w
E_w = [0 0]';
Q = [T^3/3 T^2/2;
     T^2/2 T    ] * (1/40);
P_w = Q;

% Momenty šumu v
E_v = 0;
R = 10;
P_v = R;

sim_count = 10^3; % Poèet simulací
sim_length = 26; % Délka simulací
trajectories = cell(1, sim_count); % Uložené trajektorie
measurements = cell(1, sim_count); % Uložená mìøení
KFs = cell(1, sim_count); % Uložené trajektorie odhadù Kalmanova filtru
KFs_const = cell(1, sim_count); % Uložené trajektorie odhadù 
                                    % Ustáleného Kalmanova filtru
det_est = cell(1, sim_count); % Uložené trajektorie odhadù
                              % deterministickým rekonstruktorem
KFs_inovace = zeros(sim_count, sim_length); % Uložené inovaèní posloupnosti KF
KFs_const_inovace = zeros(sim_count, sim_length); % Uložené inovaèní posloupnosti ustáleného KF

for i=1:sim_count
    trajectories{i} = zeros(2, sim_length);
    measurements{i} = zeros(1, sim_length);
end

%% SIMULACE

for iteration=1:sim_count
    x = mvnrnd(E_x0, P_x0)';
    z = C*x + mvnrnd(E_v, P_v)';
    
    trajectories{iteration}(:, 1) = x;
    measurements{iteration}(1) = z;
    
    for k=2:sim_length
        x = A*x + mvnrnd(E_w, P_w)';
        z = C*x + mvnrnd(E_v, P_v)';
        
        trajectories{iteration}(:, k) = x;
        measurements{iteration}(k) = z;
    end
end
%% KALMAN FILTER

for iteration=1:sim_count
    predict_P = P_x0;
    predict_x = E_x0;
    KFs_inovace(iteration, 1) =  measurements{iteration}(1) - predict_x(1);
    correct_P = predict_P;
    correct_x = predict_x;
    KFs{iteration}(:, 1) = correct_x;
    
    for k=2:sim_length
        predict_P = A*correct_P*A' + G*Q*G';
        predict_x = A*correct_x; % + B*u
        
        z = measurements{iteration}(k);
        KFs_inovace(iteration, k) = z - predict_x(1);
        K = predict_P*C' * ( C*predict_P*C' + R )^(-1);
        correct_P = ( eye(2) - K*C ) * predict_P * ( eye(2) - K*C )' + K*R*K';
        correct_x = predict_x + K * ( z - C*predict_x );
        
        KFs{iteration}(:, k) = correct_x;
    end
end

%% USTÁLENÁ HODNOTA KALMANOVA ZISKU

correct_P = P_x0;
K = Inf;
K_prev = -Inf;

while K_prev ~= K
    K_prev = K;
    predict_P = A*correct_P*A' + G*Q*G';
    K = predict_P*C' * ( C*predict_P*C' + R )^(-1);
    correct_P = ( eye(2) - K*C ) * predict_P * ( eye(2) - K*C )' + K*R*K';
end

K_final = K

%% USTÁLENÝ KALMANÙV FILTR
K = K_final;
for iteration=1:sim_count
    predict_P = P_x0;
    predict_x = E_x0;
    KFs_const_inovace(iteration, 1) =  measurements{iteration}(1) - predict_x(1);
    correct_P = predict_P;
    correct_x = predict_x;
    KFs_const{iteration}(:, 1) = correct_x;
    
    for k=2:sim_length
        predict_P = A*correct_P*A' + G*Q*G';
        predict_x = A*correct_x; % + B*u
        
        z = measurements{iteration}(k);
        KFs_const_inovace(iteration, k) = z - predict_x(1);
        correct_P = ( eye(2) - K*C ) * predict_P * ( eye(2) - K*C )' + K*R*K';
        correct_x = predict_x + K * ( z - C*predict_x );
        
        KFs_const{iteration}(:, k) = correct_x;
    end
end

%% DETERMINISTICKÝ REKONSTRUKTOR STAVU

K_det = [1 1]';
for iteration=1:sim_count
    estimate_x = E_x0;
    
    det_est{iteration}(:, 1) = estimate_x;
    
    for k=2:sim_length
        z = measurements{iteration}(k);
        estimate_x = (eye(2) - K_det*C)*A*estimate_x + K_det*z;
        
        det_est{iteration}(:, k) = estimate_x;
    end
end

%% POROVNÁNÍ ESTIMÁTORÙ

plot_ii(trajectories{1}, measurements{1}, KFs{1}, KFs_const{1}, det_est{1}, KFs_inovace(1, :), KFs_const_inovace(1, :)); 

%% ODHADY VARIANCÍ INOVACÍ

ks = [1 2 6 26]; % Vybrané èasové okamžiky (indexace +1)

esimate_var_KF_inovace = zeros(1, length(ks));
esimate_var_UKF_inovace = zeros(1, length(ks));
i = 1;
for k=ks
    esimate_var_KF_inovace(i) = scalar_cov(KFs_inovace(:, k), KFs_inovace(:, k));
    esimate_var_UKF_inovace(i) = scalar_cov(KFs_const_inovace(:, k), KFs_const_inovace(:, k));
    i = i + 1;
end

%% ODHADY ÈASOVÝCH KOVARIANCÍ INOVACÍ

esimate_cov_KF_inovace = zeros(length(ks), length(ks));
esimate_cov_UKF_inovace = zeros(1, length(ks));

i = 1;
for n=ks
    j = 1;
    for k=ks
        esimate_cov_KF_inovace(i, j) = scalar_cov(KFs_inovace(:, n), KFs_inovace(:, k));
        esimate_cov_UKF_inovace(i, j) = scalar_cov(KFs_const_inovace(:, n), KFs_const_inovace(:, k));
        
        j = j + 1;
    end
    i = i + 1;
end
esimate_cov_KF_inovace
esimate_cov_UKF_inovace

%% POROVNÁNÍ ODHADÙ VARIANCÍ INOVACÍ

plot_iii(esimate_var_KF_inovace, esimate_var_UKF_inovace, ks);

%% STØEDNÍCH KVADRATICKÉ CHYBY
%%% - ODHAD
e_KF = cell(1, sim_count); 
e_UKF = cell(1, sim_count); 
e_det_est = cell(1, sim_count); 

for i=1:sim_count
    e_KF{i} = ms_errors(trajectories{i}, KFs{i});
    e_UKF{i} = ms_errors(trajectories{i}, KFs_const{i});
    e_det_est{i} = ms_errors(trajectories{i}, det_est{i});
end

mean_e_KF_sim = mean_ms_errors(e_KF); % Odhad støední kvadratické chyby Kalmanova filtru
mean_e_UKF_sim = mean_ms_errors(e_UKF); % Odhad støední kvadratické chyby 
                                                          % ustáleného Kalmanova filtru
mean_e_det_est_sim = mean_ms_errors(e_det_est); % Odhad støední kvadratické chyby 
                                                % deterministického rekonstruktoru

%%% - TEORETICKÉ HODNOTY
mean_e_KF_t = cell(1, sim_length); % Teoretická hodnota støední kvadratické 
                                        % chyby Kalmanova filtru
mean_e_UKF_t = cell(1, sim_length); % Teoretická hodnota støední kvadratické chyby 
                                              % ustáleného Kalmanova filtru
mean_e_det_est_t = cell(1, sim_length); % Teoretická hodnota støední kvadratické chyby 
                                         % deterministického rekonstruktoru

correct_P = P_x0;
mean_e_KF_t{1} = correct_P;
for k=2:sim_length
    correct_P = inv( inv(A*correct_P*A' + Q) + C'*inv(R)*C );
    mean_e_KF_t{k} = correct_P;
end

correct_P = P_x0;
K = K_final;
mean_e_UKF_t{1} = correct_P;
for k=2:sim_length
    correct_P = ( eye(2) - K*C ) * (A*correct_P*A' + Q') * ( eye(2) - K*C )' + K*R*K';
    mean_e_UKF_t{k} = correct_P;
end


mean_e_det_est_t{1} = P_x0;
for k=2:sim_length
    mean_e_det_est_t{k} = (A-K_det*C*A) * mean_e_det_est_t{k-1} * (A' - A'*C'*K_det') + (eye(2) - K_det*C)*Q*(eye(2) - C'*K_det') + K_det*R*K_det';
end

plot_iv(mean_e_KF_t,mean_e_KF_sim, mean_e_UKF_t, mean_e_UKF_sim, mean_e_det_est_t, mean_e_det_est_sim)


