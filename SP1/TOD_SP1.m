clc
close all
clear all 

%% Semestralni prace c. 1 - TOD
% Jan Burian

%% I. Teorericka cast 
% Zadany system 
%x_{k+1} = F*x_k + G*w_k, w_k ~ N(0,Q) 
%z_{k} = H*x_k + v_k, v_k ~ N(0,R)

%% Ukol (i)
% Parametry systemu (zatim obecne)
syms T q R real
F = [1 T; 
     0 1]; 
G = 1; 
H = [1 0]; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 

%% Dosazeni do ziskane kovariancni matice Sigma
sigma = [R      0               0                             0; 
         0      H*Q*H' + R      H*F*Q*H'                      H*F^2*Q*H'; 
         0      H*F*Q*H'        H*F*Q*(H*F)' + H*Q*H' + R     H*F*Q*(H*F^2)' + H*F*Q*H'; 
         0      H*F^2*Q*H'      H*F*Q*(H*F^2)' + H*F*Q*H'     H*F^2*Q*(H*F^2)' + H*F*Q*(H*F)' + H*Q*H' + R];
     
%% Metoda vazenych nejmensich ctvercu
J = [H; H*F; H*F^2; H*F^3]; 
J_T = J'; 
sigma_inv = inv(sigma); 

x_mvnc = simplify(inv(J_T * sigma_inv * J) * J_T * sigma_inv);

%% Ukol(ii)
A_sigma = sigma([1,2],[1,2]); % mereni z0, z1
B_sigma = sigma([1,3],[1,3]);  % mereni z0, z2
C_sigma = sigma([2,3],[2,3]); % mereni z1, z2
D_sigma = sigma((1:3),(1:3)); % mereni z0, z1, z2

J_A = [H; H*F];
J_B = [H; H*F^2];
J_C = [H*F; H*F^2];
J_D = [H; H*F; H*F^2]; 

J_T_A = J_A'; 
J_T_B = J_B'; 
J_T_C = J_C'; 
J_T_D = J_D'; 

A_sigma_inv = inv(A_sigma); 
B_sigma_inv = inv(B_sigma); 
C_sigma_inv = inv(C_sigma); 
D_sigma_inv = inv(D_sigma);

x_mvnc_A = simplify(inv(J_T_A * A_sigma_inv * J_A) * J_T_A * A_sigma_inv); 
x_mvnc_B = simplify(inv(J_T_B * B_sigma_inv * J_B) * J_T_B * B_sigma_inv); 
x_mvnc_C = simplify(inv(J_T_C * C_sigma_inv * J_C) * J_T_C * C_sigma_inv); 
x_mvnc_D = simplify(inv(J_T_D * D_sigma_inv * J_D) * J_T_D * D_sigma_inv);

%% Dosazeni T = 1, q = 0.1, R = 1
subs_nezname = [T q R];
subs_nezname_hodnoty = [1 0.1 1];

A_odhad_subs = subs(x_mvnc_A, subs_nezname, subs_nezname_hodnoty); 
B_odhad_subs = subs(x_mvnc_B, subs_nezname, subs_nezname_hodnoty); 
C_odhad_subs = subs(x_mvnc_C, subs_nezname, subs_nezname_hodnoty); 
D_odhad_subs = subs(x_mvnc_D, subs_nezname, subs_nezname_hodnoty); 
E_odhad_subs = subs(x_mvnc, subs_nezname, subs_nezname_hodnoty);

A_sigma_subs = subs(A_sigma, subs_nezname, subs_nezname_hodnoty); 
B_sigma_subs = subs(B_sigma, subs_nezname, subs_nezname_hodnoty); 
C_sigma_subs = subs(C_sigma, subs_nezname, subs_nezname_hodnoty); 
D_sigma_subs = subs(D_sigma, subs_nezname, subs_nezname_hodnoty); 
E_sigma_subs = subs(sigma, subs_nezname, subs_nezname_hodnoty); 

%% Vypocet kovariancnich matic chyby odhadu - potrebne na vykresleni elips
% (J^{T}WJ)^{-1}

covA = inv(J_T_A * A_sigma_inv * J_A); 
covB = inv(J_T_B * B_sigma_inv * J_B);
covC = inv(J_T_C * C_sigma_inv * J_C);
covD = inv(J_T_D * D_sigma_inv * J_D);
covE = inv(J_T * sigma_inv * J);

%% Dosazeni T = 1, q = 0.1, R = 1 do matic covA, covB, covC, covD, covE
subs_nezname = [T q R];
subs_nezname_hodnoty = [1 0.1 1];

cov_A_subs = subs(covA, subs_nezname, subs_nezname_hodnoty); 
cov_B_subs = subs(covB, subs_nezname, subs_nezname_hodnoty); 
cov_C_subs = subs(covC, subs_nezname, subs_nezname_hodnoty); 
cov_D_subs = subs(covD, subs_nezname, subs_nezname_hodnoty); 
cov_E_subs = subs(covE, subs_nezname, subs_nezname_hodnoty);

%% Ukol (iii)
% Vykresleni elips
x0 = [0; 0]; 
t = 0:0.01:2*pi;

x = sin(t); 
y = cos(t); 

vector_xy = [x; y]; 
xy = zeros(2, length(t));

cells_cov = cell(1,5); 
cells_cov{1} = cov_A_subs; 
cells_cov{2} = cov_B_subs; 
cells_cov{3} = cov_C_subs; 
cells_cov{4} = cov_D_subs; 
cells_cov{5} = cov_E_subs; 

figure; 
for i = 1:length(cells_cov)
    S = cells_cov{i}; 
    P = chol(S, 'lower');
    for j = 1:length(t)
        xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
    end
    plot(xy(1,:),xy(2,:)); 

    hold on; 
end

title('3-\sigma elipsy');
xlabel('Poloha'); 
ylabel('Rychlost'); 
legend('z_0, z_1','z_0, z_2', 'z_1, z_2', 'z_0, z_1, z_2', 'z_0, z_1, z_2, z_3')
grid on; 

%% II. Simulacni cast 
%% Ukol (i)
% Vygenerovani 1000 realizaci nahodneho vektoru z 
N = 1000; 
pocet_mereni = 4; 

Z = randn(N, pocet_mereni); 

%% Vypoctene estimatory 
est_z0_z1 = double(A_odhad_subs);        
est_z0_z2 = double(B_odhad_subs);         
est_z1_z2 = double(C_odhad_subs);         
est_z0_z1_z2 = double(D_odhad_subs); 
est_z0_z1_z2_z3 = double(E_odhad_subs); 
            
%% Realizace odhadu
res_z0_z1 = zeros(2,N); %[poloha; rychlost]
res_z0_z2 = zeros(2,N); %[poloha; rychlost]
res_z1_z2 = zeros(2,N); %[poloha; rychlost]
res_z0_z1_z2 = zeros(2,N); %[poloha; rychlost]
res_z0_z1_z2_z3 = zeros(2,N); %[poloha; rychlost]

for m = 1:1:N
    mereni = Z(m,:); 
    z0 = mereni(1); 
    z1 = mereni(2);
    z2 = mereni(3);
    z3 = mereni(4);
    
    real_z0_z1 = est_z0_z1 * [z0; z1]; 
    real_z0_z2 = est_z0_z2 * [z0; z2]; 
    real_z1_z2 = est_z1_z2 * [z1; z2]; 
    real_z0_z1_z2 = est_z0_z1_z2 * [z0; z1; z2]; 
    real_z0_z1_z2_z3 = est_z0_z1_z2_z3 * [z0; z1; z2; z3]; 
    
    res_z0_z1(:,m) = real_z0_z1; 
    res_z0_z2(:,m) = real_z0_z2; 
    res_z1_z2(:,m) = real_z1_z2; 
    res_z0_z1_z2(:,m) = real_z0_z1_z2;
    res_z0_z1_z2_z3(:,m) = real_z0_z1_z2_z3;
end

%% Vykresleni realizaci 
% % z0, z1
% figure
% P = chol(cov_A_subs, 'lower');
% for j = 1:length(t)
%     xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
% end
% plot(xy(1,:),xy(2,:),'r'); 
% hold on; 
% scatter(res_z0_z1(1,:), res_z0_z1(2,:), '.', 'b');
% title('3-\sigma elipsa a realizace odhadu pro mereni z_0 a z_1');
% xlabel('Poloha'); 
% ylabel('Rychlost'); 
% grid on; 
% 
% % z0, z2
% figure
% P = chol(cov_B_subs, 'lower');
% for j = 1:length(t)
%     xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
% end
% plot(xy(1,:),xy(2,:),'r'); 
% hold on; 
% scatter(res_z0_z2(1,:), res_z0_z2(2,:), '.', 'b'); 
% title('3-\sigma elipsa a realizace odhadu pro mereni z_0 a z_2');
% xlabel('Poloha'); 
% ylabel('Rychlost'); 
% grid on; 
% 
% % z1, z2
% figure
% P = chol(cov_C_subs, 'lower');
% for j = 1:length(t)
%     xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
% end
% plot(xy(1,:),xy(2,:),'r'); 
% hold on; 
% scatter(res_z1_z2(1,:), res_z1_z2(2,:), '.', 'b'); 
% title('3-\sigma elipsa a realizace odhadu pro mereni z_1 a z_2');
% xlabel('Poloha'); 
% ylabel('Rychlost'); 
% grid on; 
% 
% % z0, z1, z2
% figure
% P = chol(cov_D_subs, 'lower');
% for j = 1:length(t)
%     xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
% end
% plot(xy(1,:),xy(2,:),'r'); 
% hold on; 
% scatter(res_z0_z1_z2(1,:), res_z0_z1_z2(2,:), '.', 'b'); 
% title('3-\sigma elipsa a realizace odhadu pro mereni z_0, z_1 a z_2');
% xlabel('Poloha'); 
% ylabel('Rychlost'); 
% grid on; 
% 
% % z0, z1, z2, z3
% figure
% P = chol(cov_E_subs, 'lower');
% for j = 1:length(t)
%     xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
% end
% plot(xy(1,:),xy(2,:),'r'); 
% hold on; 
% scatter(res_z0_z1_z2(1,:), res_z0_z1_z2(2,:), '.', 'b');
% title('3-\sigma elipsa a realizace odhadu pro mereni z_0, z_1, z_2 a z_3');
% xlabel('Poloha'); 
% ylabel('Rychlost'); 
% grid on; 

%% Histogramy - Ukol (ii)
subplot(3,2,1)
histogram(res_z0_z1(1,:), 'Normalization','pdf')
xlabel('Poloha');
title('z_0, z_1');

subplot(3,2,2)
histogram(res_z0_z2(1,:), 'Normalization','pdf');
xlabel('Poloha'); 
title('z_0, z_2');

subplot(3,2,3)
histogram(res_z1_z2(1,:), 'Normalization','pdf');
xlabel('Poloha'); 
title('z_1, z_2');

subplot(3,2,4)
histogram(res_z0_z1_z2(1,:), 'Normalization','pdf');
xlabel('Poloha'); 
title('z_0, z_1, z_2');

subplot(3,2,5)
histogram(res_z0_z1_z2_z3(1,:), 'Normalization','pdf');
xlabel('Poloha'); 
title('z_0, z_1, z_2, z_3');

%% Gaussovy krivky pro jednotliva mereni 
var_z0_z1 = cov_A_subs(1, 1); 
var_z0_z2 = cov_B_subs(1, 1);
var_z1_z2 = cov_C_subs(1, 1);
var_z0_z1_z2 = cov_D_subs(1, 1);
var_z0_z1_z2_z3 = cov_E_subs(1, 1);

mean_teo = 0; 

x_z0_z1 = sort(res_z0_z1(1,:)); 
x_z0_z2 = sort(res_z0_z2(1,:)); 
x_z1_z2 = sort(res_z1_z2(1,:)); 
x_z0_z1_z2 = sort(res_z0_z1_z2(1,:)); 
x_z0_z1_z2_z3 = sort(res_z0_z1_z2_z3(1,:)); 

y_z0_z1 = normpdf(x_z0_z1, mean_teo, sqrt(var_z0_z1));
y_z0_z2 = normpdf(x_z0_z2, mean_teo, sqrt(var_z0_z2));
y_z1_z2 = normpdf(x_z1_z2, mean_teo, sqrt(var_z1_z2));
y_z0_z1_z2 = normpdf(x_z0_z1_z2, mean_teo, sqrt(var_z0_z1_z2));
y_z0_z1_z2_z3 = normpdf(x_z0_z1_z2_z3, mean_teo, sqrt(var_z0_z1_z2_z3));

% Vykresleni 
figure; 
plot(x_z0_z1, y_z0_z1, 'LineWidth', 2); 
hold on;
histogram(res_z0_z1(1,:), 'Normalization','pdf');
title('z_0, z_1');
xlabel('Poloha'); 

figure;
plot(x_z0_z2, y_z0_z2, 'LineWidth', 2); 
hold on;
histogram(res_z0_z2(1,:), 'Normalization','pdf');
title('z_0, z_2');
xlabel('Poloha'); 

figure; 
plot(x_z1_z2, y_z1_z2, 'LineWidth', 2); 
hold on;
histogram(res_z1_z2(1,:), 'Normalization','pdf');
title('z_1, z_2');
xlabel('Poloha'); 

figure; 
plot(x_z0_z1_z2, y_z0_z1_z2, 'LineWidth', 2); 
hold on;
histogram(res_z0_z1_z2(1,:), 'Normalization','pdf');
title('z_0, z_1, z_2');
xlabel('Poloha'); 

figure; 
plot(x_z0_z1_z2_z3, y_z0_z1_z2_z3, 'LineWidth', 2); 
hold on;
histogram(res_z0_z1_z2_z3(1,:), 'Normalization','pdf');
title('z_0, z_1, z_2, z_3');
xlabel('Poloha'); 

%% Ukol (ii)
vyber_cov = cov(res_z0_z1_z2_z3'); 

figure
P = chol(cov_E_subs, 'lower');
for j = 1:length(t)
    xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
end
plot(xy(1,:),xy(2,:),'r'); 
hold on; 
P = chol(vyber_cov, 'lower');
for j = 1:length(t)
    xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
end
plot(xy(1,:),xy(2,:),'g'); 
hold on
scatter(res_z0_z1_z2(1,:), res_z0_z1_z2(2,:), '.', 'b');

title('3-\sigma elipsy a realizace odhadu pro mereni z_0, z_1, z_2 a z_3');
legend('3-\sigma elipsa (teoreticky odvozena)','3-\sigma elipsa (vyberova kovariancni matice)'); 
xlabel('Poloha'); 
ylabel('Rychlost'); 
grid on; 

