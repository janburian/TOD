clc
clear all
close all 

%% Zadany system 
%x_{k+1} = F*x_k + G*w_k, w_k ~ N(0,Q) 
%z_{k} = H*x_k + v_k, v_k ~ N(0,R)

syms T q R real
F = [1 T; 
     0 1]; 
G = 1; 
H = [1 0]; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 
      
%% Vypocet soucinu e*e^T
syms F H w0 w1 w2 w3 v0 v1 v2 v3 real

e_1  = [v0; H*w0+v1; H*F*w0+H*w1+v2; H*F^2*w0+H*F*w1+H*w2+v3]; 
chyba_1 = e_1*e_1'; 

e_2 = [v0; H*w0+v1];
chyba_2 = e_2 * e_2';

e_3 = [v0; H*F*w0+H*w1+v2];
chyba_3 = e_3 * e_3';

e_4 = [H*w0 + v1; H*F*w0 + H*w1 + v2];
chyba_4 = e_4 * e_4';

e_5 = [v0; H*w0 + v1; H*F*w0 + H*w1 + v2]
chyba_5 = e_5 * e_5'
%% Dosazeni do kovariancni matice Sigma
F = [1 T; 
     0 1]; 
G = 1; 
H = [1 0]; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 
     
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

%% Vykresleni kruznice
% t = 0:0.001:2*pi;
% 
% x = sin(t); 
% y = cos(t); 
% 
% figure; 
% plot(x,y); 

%% Vypocet kovariancnich matic chyby odhadu 
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

%% Choleskeho faktorizace
P = cov_E_subs;
S = chol(P, 'lower'); 
%chol = S*S';

%% Vykresleni elips
x0 = [0; 0]; 
t = 0:0.01:2*pi;

x = sin(t); 
y = cos(t); 

vector_xy = [x; y]; 
xy = [];

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
%%



