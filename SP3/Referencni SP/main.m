%% 3. SEMESTR�LN� PR�CE Z P�EDM�TU KKY/TOD
% autor: Tom� Honz�k
% datum: 15.12.2019

%% ZAD�N�
close all; clear; clc;

T = 1; % �asov� konstanta

A = [1 T; % Matice dynamiky
     0 1];

C = [1 0]; % Matice v�stupu

G = eye(2); % Matice �umu dynamiky

% Momenty po��te�n� podm�nky
E_x0 = [10 1]';
P_x0 = [10 0;
        0  1];

% Momenty �umu w
E_w = [0 0]';
Q = [T^3/3 T^2/2;
     T^2/2 T    ] * (1/40);
P_w = Q;

% Momenty �umu v
E_v = 0;
R = 10;
P_v = R;

sim_count = 10^3; % Po�et simulac�
sim_length = 26; % D�lka simulac�
trajectories = cell(1, sim_count); % Ulo�en� trajektorie
measurements = cell(1, sim_count); % Ulo�en� m��en�
KFs = cell(1, sim_count); % Ulo�en� trajektorie odhad� Kalmanova filtru
KFs_const = cell(1, sim_count); % Ulo�en� trajektorie odhad� 
                                    % Ust�len�ho Kalmanova filtru
det_est = cell(1, sim_count); % Ulo�en� trajektorie odhad�
                              % deterministick�m rekonstruktorem
KFs_inovace = zeros(sim_count, sim_length); % Ulo�en� inova�n� posloupnosti KF
KFs_const_inovace = zeros(sim_count, sim_length); % Ulo�en� inova�n� posloupnosti ust�len�ho KF

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

%% UST�LEN� HODNOTA KALMANOVA ZISKU

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

%% UST�LEN� KALMAN�V FILTR
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

%% DETERMINISTICK� REKONSTRUKTOR STAVU

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

%% POROVN�N� ESTIM�TOR�

plot_ii(trajectories{1}, measurements{1}, KFs{1}, KFs_const{1}, det_est{1}, KFs_inovace(1, :), KFs_const_inovace(1, :)); 

%% ODHADY VARIANC� INOVAC�

ks = [1 2 6 26]; % Vybran� �asov� okam�iky (indexace +1)

esimate_var_KF_inovace = zeros(1, length(ks));
esimate_var_UKF_inovace = zeros(1, length(ks));
i = 1;
for k=ks
    esimate_var_KF_inovace(i) = scalar_cov(KFs_inovace(:, k), KFs_inovace(:, k));
    esimate_var_UKF_inovace(i) = scalar_cov(KFs_const_inovace(:, k), KFs_const_inovace(:, k));
    i = i + 1;
end

%% ODHADY �ASOV�CH KOVARIANC� INOVAC�

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

%% POROVN�N� ODHAD� VARIANC� INOVAC�

plot_iii(esimate_var_KF_inovace, esimate_var_UKF_inovace, ks);

%% ST�EDN�CH KVADRATICK� CHYBY
%%% - ODHAD
e_KF = cell(1, sim_count); 
e_UKF = cell(1, sim_count); 
e_det_est = cell(1, sim_count); 

for i=1:sim_count
    e_KF{i} = ms_errors(trajectories{i}, KFs{i});
    e_UKF{i} = ms_errors(trajectories{i}, KFs_const{i});
    e_det_est{i} = ms_errors(trajectories{i}, det_est{i});
end

mean_e_KF_sim = mean_ms_errors(e_KF); % Odhad st�edn� kvadratick� chyby Kalmanova filtru
mean_e_UKF_sim = mean_ms_errors(e_UKF); % Odhad st�edn� kvadratick� chyby 
                                                          % ust�len�ho Kalmanova filtru
mean_e_det_est_sim = mean_ms_errors(e_det_est); % Odhad st�edn� kvadratick� chyby 
                                                % deterministick�ho rekonstruktoru

%%% - TEORETICK� HODNOTY
mean_e_KF_t = cell(1, sim_length); % Teoretick� hodnota st�edn� kvadratick� 
                                        % chyby Kalmanova filtru
mean_e_UKF_t = cell(1, sim_length); % Teoretick� hodnota st�edn� kvadratick� chyby 
                                              % ust�len�ho Kalmanova filtru
mean_e_det_est_t = cell(1, sim_length); % Teoretick� hodnota st�edn� kvadratick� chyby 
                                         % deterministick�ho rekonstruktoru

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


