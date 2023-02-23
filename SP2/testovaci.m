clc
close all
clear all 

%% Semestralni prace c. 2 - TOD
% Jan Burian

%% I. Teorericka cast 
% Zadany system 
%x_{k+1} = F*x_k + G*w_k, w_k ~ N(0,Q) 
%z_{k} = H*x_k + v_k, v_k ~ N(0,R)

%% Ukol (i)
% Parametry systemu
syms R real
q = 0.1; 
G = 1;
T = 1; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T];

%% Navrzeni estimatoru se spatnymi parametry 
% R = 1/4
clear R
R_wrong = 1/8; 
syms x_poloha x_rychlost z0 z1   
x_odhad = [x_poloha; x_rychlost]; 
H_z = [1 0; 1 1]; 
R_z_wrong_1 = [R_wrong 0; 0 R_wrong + q*((T^3)/3)]; 
Px = [1 1; 1 4]; 

P_x_z = Px * H_z'; 
P_z_wrong_R_1 = H_z * Px * H_z' + R_z_wrong_1; 
E_z = H_z * x_odhad; 
P_z_x = H_z*Px; 

LMSE_odhad_z_wrong_R_1 = x_odhad + P_x_z * inv(P_z_wrong_R_1) * ([z0; z1] - E_z); 


%% Miry duvery (kovariancni matice chyb odhadu) - nespravne parametry
format rational
% R = 1/4
mira_duvery_wrong_R_1 = Px - P_x_z * inv(P_z_wrong_R_1) * P_z_x; 

%% II. Simulacni cast
% Parametry
T = 1; 
Px = [1 1; 
      1 4];
q = 0.1; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 
R = 1; 
x0 = [0; 0]; % pocatecni podminka - radek = vektor x, sloupec =k, hloubka = iterace

F = [1 T;
     0 1];
H = [1 0]; 

%% Definice promennych 
pocet_simulaci = 1000;
Z = zeros(2, pocet_simulaci); % radek = k, sloupec = iterace

% Odhady
odhad_LMSE_z_wrong_R_1 = zeros(2,pocet_simulaci);
chyba_odhadu_LMSE_z_wrong_R_1 = zeros(2,pocet_simulaci);

%% Generovani simulaci vektoru
for i = 1:pocet_simulaci
    X(:,1,i) = x0 + (randn(1,2) * chol(Px))'; % x0
    Z(1,i) = H * X(:,1,i) + randn * sqrt(R); % z0 = Hx0 + v0
    
    X(:,2,i) = F * X(:,1,i) + (randn(1,2) * chol(Q))'; % x1 = Fx0 + w0
    Z(2,i) = H * X(:,2,i) + randn * sqrt(R); % z1 = Hx1 + v1

    % Odhady
    odhad_LMSE_z_wrong_R_1(:,i) = 1/1225 * [317*x0(1) - 120*x0(2) + 788*Z(1,i) + 120*Z(2,i); 205*x0(2) - 368*x0(1) - 652*Z(1,i) + 1020*Z(2,i)];
    
    % Vypocet chyb odhadu
    chyba_odhadu_LMSE_z_wrong_R_1(:,i) = X(:,1,i) - odhad_LMSE_z_wrong_R_1(:,i);
end

chyby_odhadu = cell(1,1); 
chyby_odhadu{1} = chyba_odhadu_LMSE_z_wrong_R_1; 


%% Histogramy 
figure; 
histogram(chyba_odhadu_LMSE_z_wrong_R_1(1,:), 'Normalization','pdf');
xlabel('Chyba odhadu polohy');
ylabel('Odhad hustoty');

%% Teoreticke Gaussovy krivky pro jednotlive pripady
var_LMSE_z_wrong_R_1 = mira_duvery_wrong_R_1(1,1);

figure; 
y_teo = normpdf(sort(chyba_odhadu_LMSE_z_wrong_R_1(1,:)), 0, sqrt(var_LMSE_z_wrong_R_1));
plot(sort(chyba_odhadu_LMSE_z_wrong_R_1(1,:)), y_teo); 
title("Teoreticky ziskane Gaussovy krivky")

%% Gaussovy krivky ze simulace pro jednotlive pripady
cell_y_sim = cell(1,1); 
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    y_sim = normpdf(sort(chyba_odhadu(1,:)), mean(chyba_odhadu(1,:)), sqrt(var(chyba_odhadu(1,:))));
    cell_y_sim{i} = y_sim; 
end

%% Variance ze simulace
for u = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{u}; 
    variance_test = var(chyba_odhadu(1,:))
end

%% Celkove porovnani krivek 
figure; 
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    p = plot(sort(chyba_odhadu(1,:)), y_teo); 
    c = p.Color; 
    hold on
    plot(sort(chyba_odhadu(1,:)), cell_y_sim{i}, '--', 'Color', c); 
    hold on
end
legend('Teoreticka', 'Simulacni')
title(sprintf('Pro falesne R = %.3f', R_wrong))



