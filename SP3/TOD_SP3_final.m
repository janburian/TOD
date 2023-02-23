clc
close all
clear all 

%% Semestralni prace c. 3 - TOD
% Jan Burian

% Zadany system 
%x_{k+1} = F*x_k + G*w_k, w_k ~ N(0,Q) 
%z_{k} = H*x_k + v_k, v_k ~ N(0,R)

%% Ukol (i)
% Parametry systemu
q = 1/40; 
T = 1;
F = [1 T; 
     0 1]; 
G = eye(2);
H = [1 0]; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 
R = 10; 

x0_mean = [10; 1];
x0_cov = [10 0; 0 1];

%% Ukol (i)
% Generovani simulaci 
pocet_trajektorii = 1000; 
pocet_kroku = 25; 

trajektorie_cell = cell(1, 1000); 
mereni_cell = cell(1, 1000); 

X = zeros(2, 26); % stavy [poloha; rychlost]
Z = zeros(26, 1); % mereni

for i = 1:pocet_trajektorii
    for k = 1:pocet_kroku+1
        if k == 1 % v prvnim kroku 
            X(:,1) = x0_mean + (randn(1,2) * chol(x0_cov))'; % x0
            Z(k,1) = H * X(:,k) + randn * sqrt(R); % z0 = Hx0 + v0
        
        else
            X(:,k) = F * X(:, k-1) + (randn(1,2) * chol(Q))'; % x_{k+1} = Fx_k + w_k
            Z(k,1) = H * X(:, k) + randn * sqrt(R); % z_k = Hx_k + v_k   
        end
    end
    trajektorie_pomocna = X'; % [poloha, rychlost]
    mereni_pomocna = Z; % [mereni]
    
    trajektorie_cell{i} = trajektorie_pomocna; 
    mereni_cell{i} = mereni_pomocna; 
end 

%% Rekonstruktor
syms k1 k2 lambda
K = [k1; k2]; 
matice = lambda * eye(2) - (eye(2) - K*H)*F; 
determinant = det(matice); 

K_subs = [1; 1]; 
rek_matice = (eye(2) - K_subs * H) * F; 

vektor_stavu_rekonstrukor = zeros(pocet_kroku+1, 2); 
rekonstruktor_result = cell(1,1000); 

x0 = [10; 1]; % pocatecni stav
for i = 1:length(trajektorie_cell)
    trajektorie = trajektorie_cell{i}; 
    mereni_vektor = mereni_cell{i}; 
    
    stav = x0; % pocatecni stav
    vektor_stavu_rekonstrukor(1,1:2) = stav;
    
    for n = 2:length(trajektorie)
        mereni = mereni_vektor(n,1); 
        stav = rek_matice * stav + K_subs * mereni; 
        vektor_stavu_rekonstrukor(n,1:2) = stav'; 
    end
    rekonstruktor_result{i} = vektor_stavu_rekonstrukor; 
end

%% Kalmanuv filtr
[Kalman_result, inovace_KF] = Kalman_filter(trajektorie_cell, mereni_cell, pocet_kroku, F, G, H, Q, R, x0_cov, x0_mean); 

%% Ziskani ustalene hodnoty Kalmanova zisku
korekce_P = x0_cov;
K = Inf; % nekonecny zisk
K_previous = -Inf;

while K_previous ~= K
    K_previous = K;
    predikce_P = F * korekce_P * F' + G * Q * G';
    K = predikce_P * H' * (H * predikce_P * H' + R)^(-1);
    korekce_P = (eye(2) - K * H) * predikce_P * (eye(2) - K * H)' + K * R * K';
end

K_inf = K; 

%% Ustaleny Kalmanuv filtr
[Kalman_ust_result, inovace_KF_ust] = Kalman_filter_ust(trajektorie_cell, mereni_cell, pocet_kroku, F, G, H, Q, R, x0_cov, x0_mean, K_inf);


%% Ukol (ii)
trajektorie_1 = trajektorie_cell{1}; 
mereni_1 = mereni_cell{1}; 

rekonstruktor_trajektorie_1 = rekonstruktor_result{1};
KF_trajektorie_1 = Kalman_result{1}; 
KF_ust_trajektorie_1 = Kalman_ust_result{1}; 
 
% Porovnani odhadu polohy
figure
hold on 
plot(0:length(trajektorie_1)-1, trajektorie_1(:, 1))
plot(0:length(rekonstruktor_trajektorie_1)-1, rekonstruktor_trajektorie_1(:, 1)); 
plot(0:length(KF_trajektorie_1)-1, KF_trajektorie_1(:, 1));
plot(0:length(KF_ust_trajektorie_1)-1, KF_ust_trajektorie_1(:, 1));
mereni_1(1) = NaN; 
plot(0:length(mereni_1)-1, mereni_1(:, 1), 'o', 'LineWidth', 1.2);

xlabel('k'); 
ylabel('Poloha'); 
legend('Skutecna trajektorie', 'Rekonstruktor', 'Kalmanuv filtr', 'Ustaleny Kalmanuv filtr',  'Mereni'); 
title('Porovnani odhadu polohy');  


% Porovnani odhadu rychlosti
figure 
hold on 
plot(0:length(trajektorie_1)-1, trajektorie_1(:, 2))
plot(0:length(rekonstruktor_trajektorie_1)-1, rekonstruktor_trajektorie_1(:, 2));
plot(0:length(KF_trajektorie_1)-1, KF_trajektorie_1(:, 2));
plot(0:length(KF_ust_trajektorie_1)-1, KF_ust_trajektorie_1(:, 2));

xlabel('k'); 
ylabel('Rychlost'); 
legend('Skutecna trajektorie', 'Rekonstruktor', 'Kalmanuv filtr', 'Ustaleny Kalmanuv filtr'); 
title('Porovnani odhadu rychlosti');  

% Chyby odhadu
chyba_rekonstruktor = trajektorie_1 - rekonstruktor_trajektorie_1; 
chyba_KF = trajektorie_1 - KF_trajektorie_1;  
chyba_KF_ust = trajektorie_1 - KF_ust_trajektorie_1; 

% Chyby odhadu polohy
colorOrder = get(gca, 'ColorOrder');

figure 
hold on 
plot(0:length(rekonstruktor_trajektorie_1)-1, chyba_rekonstruktor(:, 1), 'Color', colorOrder(2,:));
plot(0:length(KF_trajektorie_1)-1, chyba_KF(:, 1), 'Color', colorOrder(3,:));
plot(0:length(KF_ust_trajektorie_1)-1, chyba_KF_ust(:, 1), 'Color', colorOrder(4,:));
plot(0:length(KF_trajektorie_1)-1, zeros(1, length(KF_trajektorie_1)), 'b--');

xlabel('k'); 
ylabel('Chyba odhadu polohy'); 
legend('Rekonstruktor', 'Kalmanuv filtr', 'Ustaleny Kalmanuv filtr'); 
title('Porovnani chyb odhadu polohy'); 

% Chyby odhadu rychlosti
figure
hold on 
plot(0:length(rekonstruktor_trajektorie_1)-1, chyba_rekonstruktor(:, 2), 'Color', colorOrder(2,:));
plot(0:length(KF_trajektorie_1)-1, chyba_KF(:, 2), 'Color', colorOrder(3,:));
plot(0:length(KF_ust_trajektorie_1)-1, chyba_KF_ust(:, 2), 'Color', colorOrder(4,:));
plot(0:length(KF_trajektorie_1)-1, zeros(1, length(KF_trajektorie_1)), 'b--');

xlabel('k'); 
ylabel('Chyba odhadu rychlosti'); 
legend('Rekonstruktor', 'Kalmanuv filtr', 'Ustaleny Kalmanuv filtr'); 
title('Porovnani chyb odhadu rychlosti'); 

%% Ukol (iii)
colorOrder = get(gca, 'ColorOrder');
% figure
% hold on 
% plot(0:length(inovace_KF(:,1))-1, inovace_KF(:, 1), 'Color', colorOrder(3,:));
% plot(0:length(inovace_KF_ust(:,1))-1, inovace_KF_ust(:, 1), 'Color', colorOrder(4,:));
% xlabel('k'); 
% ylabel('Hodnota inovace'); 
% legend('Kalmanuv filtr', 'Ustaleny Kalmanuv filtr'); 
% title('Porovnani inovacnich posloupnosti pro prvni trajektorii');

result_means_inovace_KF = zeros(pocet_kroku+1, 1); 
result_means_inovace_KF_ust = zeros(pocet_kroku+1, 1); 

for i = 1:pocet_kroku+1
    inovace_KF_ki = inovace_KF(i,:); 
    result_means_inovace_KF(i, 1) = mean(inovace_KF_ki);
    
    inovace_KF_ust_ki = inovace_KF_ust(i, :); 
    result_means_inovace_KF_ust(i, 1) = mean(inovace_KF_ust_ki); 
end

figure; 
hold on 
plot(0:length(result_means_inovace_KF)-1, result_means_inovace_KF(:, 1), 'Color', colorOrder(3,:));
plot(0:length(result_means_inovace_KF_ust)-1, result_means_inovace_KF_ust(:, 1), 'Color', colorOrder(4,:));
xlabel('k'); 
ylabel('Odhad stredni hodnoty inovace'); 
legend('Kalmanuv filtr', 'Ustaleny Kalmanuv filtr'); 
title('Porovnani odhadu strednich hodnot inovaci');

vyber_kroky = [1 2 6 26]; % k = 0, 1, 5 a 26 (posunute o jeden index dopredu)
kov_matice_inovace_KF = zeros(length(vyber_kroky), length(vyber_kroky)); 
kov_matice_inovace_KF_ust = zeros(length(vyber_kroky), length(vyber_kroky));

% Kovariance
i = 0; 
j = 0; 
for k = vyber_kroky
    clen1_KF = inovace_KF(k, :); 
    clen1_KF_ust = inovace_KF_ust(k, :); 
    
    i = mod(i, length(vyber_kroky))+1; 
    for n = vyber_kroky
        j = mod(j, length(vyber_kroky))+1; 
        
        clen2_KF = inovace_KF(n, :); 
        clen2_KF_ust = inovace_KF_ust(n, :);  
        
        kovariance_KF = cov(clen1_KF, clen2_KF); 
        kovariance_KF_ust = cov(clen1_KF_ust, clen2_KF_ust); 
        
        kov_matice_inovace_KF(i, j) = kovariance_KF(2,1); % kovariance = [Px Pxz, Pzx Pz]
        kov_matice_inovace_KF_ust(i, j) = kovariance_KF_ust(2,1); % kovariance = [Px Pxz, Pzx Pz]
    end
end

% Variance (diagonaly z kovariancnich matic)
varis_inovace_KF = diag(kov_matice_inovace_KF); 
varis_inovace_KF_ust = diag(kov_matice_inovace_KF_ust); 

figure; 
hold on
plot(vyber_kroky-1, varis_inovace_KF, 'o', 'Color', colorOrder(3,:), 'MarkerSize', 10);
plot(vyber_kroky-1, varis_inovace_KF_ust, '+', 'Color', colorOrder(4,:), 'MarkerSize', 10); 
title('Porovnani odhadu varianci inovaci'); 
xlabel('k'); 
ylabel('Odhad variance inovace')
legend('Kalmanuv filtr', 'Ustaleny Kalmanuv filtr')

%% Ukol (iv)
% Odhady strednich kvadratickych chyb (MSE) (ze simulace)
MSE_rekonstruktor_pomocna = cell(1, pocet_trajektorii); 
MSE_Kalman_pomocna = cell(1, pocet_trajektorii); 
MSE_Kalman_ust_pomocna = cell(1, pocet_trajektorii);

% Ziskani matic z cell
for i = 1:pocet_trajektorii
    trajektorie = trajektorie_cell{i}; 
    odhad_rekonstruktor = rekonstruktor_result{i}; 
    odhad_Kalman = Kalman_result{i}; 
    odhad_Kalman_ust = Kalman_ust_result{i}; 
    
    MSE_rekonstruktor_pomocna{i} = count_MSE_pomocna(trajektorie, odhad_rekonstruktor); 
    MSE_Kalman_pomocna{i} = count_MSE_pomocna(trajektorie, odhad_Kalman); 
    MSE_Kalman_ust_pomocna{i} = count_MSE_pomocna(trajektorie, odhad_Kalman_ust); 
end

MSE_rekonstruktor = count_MSE(MSE_rekonstruktor_pomocna); 
MSE_Kalman = count_MSE(MSE_Kalman_pomocna); 
MSE_Kalman_ust = count_MSE(MSE_Kalman_ust_pomocna); 

% Odhady strednich kvadratickych chyb (MSE) (teoreticke)
MSE_rekonstruktor_teo = cell(1, pocet_kroku+1); 
MSE_Kalman_teo = cell(1, pocet_kroku+1); 
MSE_Kalman_ust_teo = cell(1, pocet_kroku+1); 

MSE_rekonstruktor_teo{1} = x0_cov;
K = [1; 1]; 
for k = 2:pocet_kroku+1
    MSE_rekonstruktor_teo{k} = (F - K * H * F) * MSE_rekonstruktor_teo{k-1} * (F' - F' * H' * K') + (eye(2) - K * H)* Q * (eye(2) - H'*K') + K * R * K';
end

korekce_P = x0_cov;
MSE_Kalman_teo{1} = korekce_P;
for k = 2:pocet_kroku+1
    korekce_P = inv(inv(F * korekce_P * F' + Q) + H' * inv(R) * H);
    MSE_Kalman_teo{k} = korekce_P;
end

korekce_P = x0_cov;
K = K_inf;
MSE_Kalman_ust_teo{1} = korekce_P;
for k = 2:pocet_kroku+1
    korekce_P = (eye(2) - K * H) * (F * korekce_P * F' + Q') * (eye(2) - K * H)' + K * R * K';
    MSE_Kalman_ust_teo{k} = korekce_P;
end

% Prevod na matice z cell na matice
MSE_rekonstruktor_teo_matice = [MSE_rekonstruktor_teo{:}];
MSE_Kalman_teo_matice = [MSE_Kalman_teo{:}]; 
MSE_Kalman_ust_teo_matice = [MSE_Kalman_ust_teo{:}]; 


% Vykresleni 
colorOrder = get(gca, 'ColorOrder');

figure; 
hold on; 
plot(0:length(MSE_rekonstruktor(:,1))-1, MSE_rekonstruktor(:, 1), 'Color', colorOrder(2,:));
plot(0:length(MSE_rekonstruktor(:,1))-1, MSE_rekonstruktor_teo_matice(1,1:2:end), ':', 'Color', colorOrder(2,:)); 

plot(0:length(MSE_Kalman(:,1))-1, MSE_Kalman(:, 1), 'Color', colorOrder(3,:));
plot(0:length(MSE_Kalman(:,1))-1, MSE_Kalman_teo_matice(1, 1:2:end), ':', 'Color', colorOrder(3,:));

plot(0:length(MSE_Kalman_ust(:,1))-1, MSE_Kalman_ust(:, 1), 'Color', colorOrder(4,:));
plot(0:length(MSE_Kalman_ust(:,1))-1, MSE_Kalman_ust_teo_matice(1, 1:2:end), ':', 'Color', colorOrder(4,:));

xlabel('k'); 
ylabel('Stredni kvadraticka chyba polohy'); 
legend('Rekonstruktor - simulace', 'Rekonstruktor - teorie', 'Kalmanuv filtr - simulace', 'Kalmanuv filtr - teorie', 'Ustaleny Kalmanuv filtr - simulace', 'Ustaleny Kalmanuv filtr - teorie');  
title('Porovnani MSE (poloha)');

figure; 
hold on 
plot(0:length(MSE_rekonstruktor(:,1))-1, MSE_rekonstruktor(:, 2), 'Color', colorOrder(2,:));
plot(0:length(MSE_rekonstruktor(:,1))-1, MSE_rekonstruktor_teo_matice(2, 2:2:end), ':', 'Color', colorOrder(2,:)); 

plot(0:length(MSE_Kalman(:,1))-1, MSE_Kalman(:, 2), 'Color', colorOrder(3,:));
plot(0:length(MSE_Kalman(:,1))-1, MSE_Kalman_teo_matice(2, 2:2:end), ':', 'Color', colorOrder(3,:));

plot(0:length(MSE_Kalman_ust(:,1))-1, MSE_Kalman_ust(:, 2), 'Color', colorOrder(4,:));
plot(0:length(MSE_Kalman_ust(:,1))-1, MSE_Kalman_ust_teo_matice(2, 2:2:end), ':', 'Color', colorOrder(4,:));
xlabel('k'); 
ylabel('Stredni kvadraticka chyba rychlosti'); 
legend('Rekonstruktor - simulace', 'Rekonstruktor - teorie', 'Kalmanuv filtr - simulace', 'Kalmanuv filtr - teorie', 'Ustaleny Kalmanuv filtr - simulace', 'Ustaleny Kalmanuv filtr - teorie'); 
title('Porovnani MSE (rychlost)');

% Stredni kvadraticka chyba
% Simulace
MSE_rekonstruktor_sim_stopy = MSE_rekonstruktor(:, 1) + MSE_rekonstruktor(:, 2); 
MSE_Kalman_sim_stopy = MSE_Kalman(:, 1) + MSE_Kalman(:, 2); 
MSE_Kalman_ust_sim_stopy = MSE_Kalman_ust(:, 1) + MSE_Kalman_ust(:, 2);

% Teorie 
MSE_rekonstruktor_teo_stopy = count_traces(MSE_rekonstruktor_teo); 
MSE_Kalman_teo_stopy = count_traces(MSE_Kalman_teo);
MSE_Kalman_ust_teo_stopy = count_traces(MSE_Kalman_ust_teo); 


figure; 
hold on 
plot(0:length(MSE_rekonstruktor_sim_stopy(:,1))-1, MSE_rekonstruktor_sim_stopy(:, 1), 'Color', colorOrder(2,:));
plot(0:length(MSE_rekonstruktor_teo_stopy(:,1))-1, MSE_rekonstruktor_teo_stopy(:, 1), ':', 'Color', colorOrder(2,:)); 

plot(0:length(MSE_Kalman_sim_stopy(:,1))-1, MSE_Kalman_sim_stopy(:, 1), 'Color', colorOrder(3,:));
plot(0:length(MSE_Kalman_teo_stopy(:,1))-1, MSE_Kalman_teo_stopy(:, 1), ':', 'Color', colorOrder(3,:));

plot(0:length(MSE_Kalman_ust_sim_stopy(:,1))-1, MSE_Kalman_ust_sim_stopy(:, 1), 'Color', colorOrder(4,:));
plot(0:length(MSE_Kalman_ust_teo_stopy(:,1))-1, MSE_Kalman_ust_teo_stopy(:, 1), ':', 'Color', colorOrder(4,:));
xlabel('k'); 
ylabel('Stredni kvadraticka chyba'); 
legend('Rekonstruktor - simulace', 'Rekonstruktor - teorie', 'Kalmanuv filtr - simulace', 'Kalmanuv filtr - teorie', 'Ustaleny Kalmanuv filtr - simulace', 'Ustaleny Kalmanuv filtr - teorie'); 
title('Porovnani MSE');

    
