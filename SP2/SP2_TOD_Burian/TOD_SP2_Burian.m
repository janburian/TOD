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
syms T R real
q = 0.1; 
F = [1 T; 
     0 1]; 
G = 1; 
H = [1 0]; 
Q = q * [T^3/3   T^2/2; 
         T^2/2   T]; 
R = 1; 
%% Metoda maximalni verohodnosti pro mereni z_0, z_1
H = [1 0; 1 T]; 
sigma = [R 0; 0 (q*(T^3)/ 3) + R]; 

ML_odhad = inv(H' * inv(sigma) * H) * H' * inv(sigma); 
ML_odhad_subs = subs(ML_odhad, T, 1);

% kovariancni matice (chyby) odhadu
cov_odhad_ML = inv(H' * inv(sigma) * H);
cov_odhad_ML_subs = subs(cov_odhad_ML, T, 1); 

%% Linearni odhady ve smyslu stredni kvadraticke chyby pro jednotliva mereni
% mereni z0
syms x_poloha x_rychlost z0   
x_odhad = [x_poloha; x_rychlost]; 
H_z0 = [1 0]; 
R_z0 = R; 
Px = [1 1; 1 4]; 

P_x_z0 = Px * H_z0'; 
P_z0 = H_z0 * Px * H_z0' + R_z0; 
E_z0 = H_z0 * x_odhad; 

LMSE_odhad_z0 = x_odhad + P_x_z0 * inv(P_z0) * (z0 - E_z0); 

% mereni z1
syms x_poloha x_rychlost z1   
T = 1; 

x_odhad = [x_poloha; x_rychlost]; 
H_z1 = [1 T]; 
R_z1 = R + q*((T^3)/3); 
Px = [1 1; 1 4]; 

P_x_z1 = Px * H_z1'; 
P_z1 = H_z1 * Px * H_z1' + R_z1; 
E_z1 = H_z1 * x_odhad; 

LMSE_odhad_z1 = x_odhad + P_x_z1 * inv(P_z1) * (z1 - E_z1); 

% mereni z = [z0; z1]
syms x_poloha x_rychlost z0 z1   
x_odhad = [x_poloha; x_rychlost]; 
H_z = [1 0; 1 1]; 
R_z = [R 0; 0 R + q*((T^3)/3)]; 
Px = [1 1; 1 4]; 

P_x_z = Px * H_z'; 
P_z = H_z * Px * H_z' + R_z; 
E_z = H_z * x_odhad; 

LMSE_odhad_z = x_odhad + P_x_z * inv(P_z) * ([z0; z1] - E_z); 

%% Miry duvery (kovariancni matice chyb odhadu)
format rational
% mereni z0
P_z0_x = H_z0*Px; 
mira_duvery_z0 = Px - P_x_z0 * inv(P_z0) * P_z0_x; 

% mereni z1
P_z1_x = H_z1*Px; 
mira_duvery_z1 = Px - P_x_z1 * inv(P_z1) * P_z1_x; 

% mereni z
P_z_x = H_z*Px; 
mira_duvery_z = Px - P_x_z * inv(P_z) * P_z_x; 

%% Kovariancni matice odhadu
cov_LMSE_z0 = P_x_z0 * inv(P_z0) * P_z0_x; 
cov_LMSE_z1 = P_x_z1 * inv(P_z1) * P_z1_x; 
cov_LMSE_z = P_x_z * inv(P_z) * P_z_x; 

%% Navrzeni estimatoru se spatnymi parametry 
% R = 1/4
clear R
R = 1/4; 
syms x_poloha x_rychlost z0 z1   
x_odhad = [x_poloha; x_rychlost]; 
H_z = [1 0; 1 1]; 
R_z_wrong_1 = [R 0; 0 R + q*((T^3)/3)]; 
Px = [1 1; 1 4]; 

P_x_z = Px * H_z'; 
P_z_wrong_R_1 = H_z * Px * H_z' + R_z_wrong_1; 
E_z = H_z * x_odhad; 

LMSE_odhad_z_wrong_R_1 = x_odhad + P_x_z * inv(P_z_wrong_R_1) * ([z0; z1] - E_z); 

% R = 4 
clear R
R = 4; 
syms x_poloha x_rychlost z0 z1   
x_odhad = [x_poloha; x_rychlost]; 
H_z = [1 0; 1 1]; 
R_z_wrong_2 = [R 0; 0 R + q*((T^3)/3)]; 
Px = [1 1; 1 4]; 

P_x_z = Px * H_z'; 
P_z_wrong_R_2 = H_z * Px * H_z' + R_z_wrong_2; 
E_z = H_z * x_odhad; 

LMSE_odhad_z_wrong_R_2 = x_odhad + P_x_z * inv(P_z_wrong_R_2) * ([z0; z1] - E_z); 

% x0^p := x0^p - 5
clear R
R = 1; 
x_odhad = [x_poloha-5; x_rychlost]; 
H_z = [1 0; 1 1]; 
R_z = [R 0; 0 R + q*((T^3)/3)]; 
Px = [1 1; 1 4]; 

P_x_z = Px * H_z'; 
P_z = H_z * Px * H_z' + R_z; 
E_z = H_z * x_odhad; 

LMSE_odhad_z_posun_x_poloha = x_odhad + P_x_z * inv(P_z) * ([z0; z1] - E_z); 

%% Strannost odhadu 
E_x = [x_poloha; x_rychlost]; 
E_x_odhad = [x_poloha-5; x_rychlost]; % odhad s posunutou polohou o -5

strannost = E_x - E_x_odhad - P_x_z * inv(P_z) * H_z * E_x + P_x_z * inv(P_z) * H_z * E_x_odhad; 

%% Miry duvery (kovariancni matice chyb odhadu) - nespravne parametry
format rational
% R = 1/4
mira_duvery_wrong_R_1 = Px - P_x_z * inv(P_z_wrong_R_1) * P_z_x; 

% R = 4 
mira_duvery_wrong_R_2 = Px - P_x_z * inv(P_z_wrong_R_2) * P_z_x; 

% x0^p := x0^p - 5 
mira_duvery_posun_x_poloha = Px - P_x_z * inv(P_z) * P_z_x; 

%% Vykresleni elips pro kovariancni matice chyb odhadu
x0 = [0; 0]; 
t = 0:0.01:2*pi;

x = sin(t); 
y = cos(t); 

vector_xy = [x; y]; 
xy = zeros(2, length(t));

cells_cov = cell(1,4); 
cells_cov{1} = cov_odhad_ML_subs; 
cells_cov{2} = mira_duvery_z0; 
cells_cov{3} = mira_duvery_z1; 
cells_cov{4} = mira_duvery_z; 

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
scatter(x0(1), x0(2), '+', 'k');

title('3-\sigma elipsy odhadu - ML a LMSE');
xlabel('Odhad polohy'); 
ylabel('Odhad rychlosti'); 
legend('ML: $z$','LMSE: $z_0$', 'LMSE: $z_1$', 'LMSE: $z$', 'Interpreter', 'latex')
grid on; 

%% Vykresleni elips pro kovariancni matice chyb odhadu zalozene na nespr. hodnotach parametru
cells_cov = cell(1,4); 
cells_cov{1} = mira_duvery_z; 
cells_cov{2} = mira_duvery_wrong_R_1; 
cells_cov{3} = mira_duvery_wrong_R_2; 
cells_cov{4} = mira_duvery_posun_x_poloha; 

figure; 
for i = 1:length(cells_cov)
    S = cells_cov{i}; 
    P = chol(S, 'lower');
    if i == 4
        x0 = [2.5; -2.5]; 
        for j = 1:length(t)
            xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
        end
        plot(xy(1,:), xy(2,:));
        scatter(x0(1), x0(2), '+', 'k'); 
        hold on; 
    else 
        x0 = [0; 0]; 
        for j = 1:length(t)
            xy(:,j) = x0 + 3 * P * vector_xy(:,j); 
        end
        plot(xy(1,:), xy(2,:));
        hold on;
    end
end
scatter(0, 0, '+', 'k');

title('3-\sigma elipsy odhadu LMSE (falesne a skutecne parametry)');
xlabel('Odhad polohy'); 
ylabel('Odhad rychlosti'); 
legend('LMSE: $z$','LMSE: $R = \frac{1}{4}$', 'LMSE: $R = 4$', 'LMSE: $x_0^p = x_0^p - 5$', 'Interpreter', 'latex')
grid on; 


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
odhad_ML_z = zeros(2,pocet_simulaci);

odhad_LMSE_z0 = zeros(2,pocet_simulaci);
odhad_LMSE_z1 = zeros(2,pocet_simulaci);
odhad_LMSE_z = zeros(2,pocet_simulaci);

odhad_LMSE_z_wrong_R_1 = zeros(2,pocet_simulaci);
odhad_LMSE_z_wrong_R_2 = zeros(2,pocet_simulaci);
odhad_LMSE_z_posun_x_poloha = zeros(2,pocet_simulaci);

% Chyby odhadu
chyba_odhadu_ML_z = zeros(2,pocet_simulaci);

chyba_odhadu_LMSE_z0 = zeros(2,pocet_simulaci);
chyba_odhadu_LMSE_z1 = zeros(2,pocet_simulaci);
chyba_odhadu_LMSE_z = zeros(2,pocet_simulaci);

chyba_odhadu_LMSE_z_wrong_R_1 = zeros(2,pocet_simulaci);
chyba_odhadu_LMSE_z_wrong_R_2 = zeros(2,pocet_simulaci);
chyba_odhadu_LMSE_z_posun_x_poloha = zeros(2,pocet_simulaci);

%% Generovani simulaci vektoru
for i = 1:pocet_simulaci
    X(:,1,i) = x0 + (randn(1,2) * chol(Px))'; % x0
    Z(1,i) = H * X(:,1,i) + randn * sqrt(R); % z0 = Hx0 + v0
    
    X(:,2,i) = F * X(:,1,i) + (randn(1,2) * chol(Q))'; % x1 = Fx0 + w0
    Z(2,i) = H * X(:,2,i) + randn * sqrt(R); % z1 = Hx1 + v1

    % Odhady
    odhad_ML_z(:,i) = [1 0; -1 1] * Z([1,2],i);
    
    odhad_LMSE_z0(:,i) = [x0(1)/2 + Z(1,i)/2; x0(2) - x0(1)/2 + Z(1,i)/2];
    odhad_LMSE_z1(:,i) = 1/241 * [181*x0(1) - 60*x0(2) + 60*Z(2,i); 91*x0(2) - 150*x0(1) + 150*Z(2,i)];
    odhad_LMSE_z(:,i) = 1/362 * [181*x0(1) - 60*x0(2) + 121*Z(1,i) + 60*Z(2,i); 122*x0(2) - 181*x0(1) - 59*Z(1,i) + 240*Z(2,i)];
    
    odhad_LMSE_z_wrong_R_1(:,i) = 1/1225 * [317*x0(1) - 120*x0(2) + 788*Z(1,i) + 120*Z(2,i); 205*x0(2) - 368*x0(1) - 652*Z(1,i) + 1020*Z(2,i)];
    odhad_LMSE_z_wrong_R_2(:,i) = 1/1535 * [1084*x0(1) - 240*x0(2) + 211*Z(1,i) + 240*Z(2,i); 845*x0(2) - 721*x0(1) + 31*Z(1,i) + 690*Z(2,i)];
    odhad_LMSE_z_posun_x_poloha(:,i) = 1/362 * [181*x0(1) - 60*x0(2) + 121*Z(1,i) + 60*Z(2,i) - 905; 122*x0(2) - 181*x0(1) - 59*Z(1,i) + 240*Z(2,i) + 905];
    
    % Vypocet chyb odhadu
    chyba_odhadu_ML_z(:,i) = X(:,1,i) - odhad_ML_z(:,i);
    chyba_odhadu_LMSE_z0(:,i) = X(:,1,i) - odhad_LMSE_z0(:,i);
    chyba_odhadu_LMSE_z1(:,i) = X(:,1,i) - odhad_LMSE_z1(:,i);
    chyba_odhadu_LMSE_z(:,i) = X(:,1,i) - odhad_LMSE_z(:,i);
    chyba_odhadu_LMSE_z_wrong_R_1(:,i) = X(:,1,i) - odhad_LMSE_z_wrong_R_1(:,i);
    chyba_odhadu_LMSE_z_wrong_R_2(:,i) = X(:,1,i) - odhad_LMSE_z_wrong_R_2(:,i);
    chyba_odhadu_LMSE_z_posun_x_poloha(:,i) = X(:,1,i) - odhad_LMSE_z_posun_x_poloha(:,i);
end

chyby_odhadu = cell(1,7); 
chyby_odhadu{1} = chyba_odhadu_ML_z; 
chyby_odhadu{2} = chyba_odhadu_LMSE_z0; 
chyby_odhadu{3} = chyba_odhadu_LMSE_z1; 
chyby_odhadu{4} = chyba_odhadu_LMSE_z; 
chyby_odhadu{5} = chyba_odhadu_LMSE_z_wrong_R_1; 
chyby_odhadu{6} = chyba_odhadu_LMSE_z_wrong_R_2; 
chyby_odhadu{7} = chyba_odhadu_LMSE_z_posun_x_poloha; 

%% Histogramy 
nbins = 30; 
titles = ["ML: $z$", "LMSE: $z_0$", "LMSE: $z_1$", "LMSE: $z$", "LMSE: $z; R = \frac{1}{4}$", "LMSE: $z; R = 4$", "LMSE: $z$; $x_0^{p} = x_0^{p} - 5$"]; 
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    figure; 
    histogram(chyba_odhadu(1,:), nbins, 'Normalization','pdf');
    xlabel('Chyba odhadu polohy');
    ylabel('Odhad hustoty');
    title(titles(i), 'Interpreter','latex');
end

figure; 
for j = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{j};
    subplot(4,2,j)
    histogram(chyba_odhadu(1,:), nbins, 'Normalization','pdf')
    xlabel('Chyba odhadu polohy');
    ylabel('Odhad hustoty')
    title(titles(j), 'Interpreter','latex');
    xlim([-3 5])
    ylim([0 0.8])
end

%% Teoreticke Gaussovy krivky pro jednotlive pripady
var_ML_z = cov_odhad_ML_subs(1,1); 
var_LMSE_z0 = mira_duvery_z0(1,1);
var_LMSE_z1 = mira_duvery_z1(1,1);
var_LMSE_z = mira_duvery_z(1,1);
var_LMSE_z_wrong_R_1 = mira_duvery_wrong_R_1(1,1);
var_LMSE_z_wrong_R_2 = mira_duvery_wrong_R_2(1,1);
var_LMSE_z_posun_x_poloha = mira_duvery_posun_x_poloha(1,1); 

cell_varis = cell(1,7); 
cell_varis{1} = var_ML_z; 
cell_varis{2} = var_LMSE_z0; 
cell_varis{3} = var_LMSE_z1; 
cell_varis{4} = var_LMSE_z; 
cell_varis{5} = var_LMSE_z_wrong_R_1; 
cell_varis{6} = var_LMSE_z_wrong_R_2; 
cell_varis{7} = var_LMSE_z_posun_x_poloha; 

mean_teo = [0 0 0 0 0 0 2.5]; 
cell_y_teo = cell(1,7); 
figure; 
for j = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{j}; 
    variance = cell_varis{j}; 
    y_teo = normpdf(sort(chyba_odhadu(1,:)), mean_teo(j), sqrt(variance));
    cell_y_teo{j} = y_teo; 
    plot(sort(chyba_odhadu(1,:)), y_teo); 
    hold on;
end
title("Porovnani teoretickych hustot"); 
xlabel('Chyba odhadu polohy');
ylabel('Odhad hustoty');
legend(titles, 'Interpreter', 'latex')

%% Gaussovy krivky ze simulace pro jednotlive pripady
cell_y_sim = cell(1,7); 
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    y_sim = normpdf(sort(chyba_odhadu(1,:)), mean(chyba_odhadu(1,:)), sqrt(var(chyba_odhadu(1,:))));
    cell_y_sim{i} = y_sim; 
end

%% Variance ze simulace
for u = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{u}; 
    variance_test = var(chyba_odhadu(1,:)); 
end

%% Porovnani jednotlivych krivek mezi sebou 
for j = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{j}; 
    figure; 
    plot(sort(chyba_odhadu(1,:)), cell_y_teo{j}); 
    hold on;
    plot(sort(chyba_odhadu(1,:)), cell_y_sim{j}); 
end

%% Porovnani prolozenych Gaussovych krivek mezi sebou 
figure; 
for j = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{j};  
    plot(sort(chyba_odhadu(1,:)), cell_y_sim{j}); 
    hold on;
end
title("Porovnani prolozenych hustot");
xlabel('Chyba odhadu polohy');
ylabel('Odhad hustoty');
legend(titles, 'Interpreter', 'latex')
ylim([0 1])

%% Gaussovy krivky (metoda prirazeni momentu) a histogramy 
%titles = ["ML (mereni z)", "LMSE (mereni z_0)", "LMSE (mereni z_1)", "LMSE (mereni z)", "LMSE (mereni z s parametrem R = 1/4)", "LMSE (mereni z s parametrem R = 4)", "LMSE (mereni z s posunutym parametrem x_0^{poloha})"]; 
for i = 1:length(chyby_odhadu)
    figure; 
    chyba_odhadu = chyby_odhadu{i};  
    y_sim = cell_y_sim{i}; 
    plot(sort(chyba_odhadu(1,:)), y_sim, 'LineWidth', 1.5); 
    hold on;
    histogram(chyba_odhadu(1,:), nbins, 'Normalization','pdf');
    xlabel('Chyba odhadu polohy');
    ylabel('Odhad hustoty');
    title(titles(i), 'Interpreter', 'latex');
end

%% Celkove porovnani krivek (prolozene a teoreticke)
figure; 
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    p = plot(sort(chyba_odhadu(1,:)), cell_y_teo{i}); 
    c = p.Color; 
    hold on
    plot(sort(chyba_odhadu(1,:)), cell_y_sim{i}, '--', 'Color', c); 
    hold on
end

legend('ML: $z$ (teor)', 'ML: $z$ (sim)', 'LMSE: $z_0$ (teor)', 'LMSE: $z_0$ (sim)', 'LMSE: $z_1$ (teor)', 'LMSE: $z_1$ (sim)', 'LMSE: $z$ (teor)', 'LMSE: $z$ (sim)', 'LMSE: $z; R = \frac{1}{4}$ (teor)', 'LMSE: $z; R = \frac{1}{4}$ (sim)', 'LMSE: $z; R = 4$ (teor)', 'LMSE: $z; R = 4$ (sim)', 'LMSE: $z$; posun $x_0^{p}$ (teor)', 'LMSE: $z$; posun $x_0^{p}$ (sim)', 'Interpreter', 'latex')
title("Porovnani prolozenych a teoretickych hustot"); 
xlabel('Chyba odhadu polohy');
ylabel('Odhad hustoty');

%% Histogramy 
figure
for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    histogram(chyba_odhadu(1,:), nbins, 'Normalization','pdf', 'DisplayStyle', 'stairs');
    %histogram(chyba_odhadu(1,:), nbins, 'Normalization','pdf');
    hold on
    xlabel('Chyba odhadu polohy');
    ylabel('Odhad hustoty');
end
legend("ML: $z$", "LMSE: $z_0$", "LMSE: $z_1$", "LMSE: $z$", "LMSE: $z; R = \frac{1}{4}$", "LMSE: $z; R = 4$", "LMSE: $z$; $x_0^{p} = x_0^{p} - 5$", 'Location','northwest', 'Interpreter','latex'); 
title("Porovnani jednotlivych histogramu");

%% Stredni hodnoty a variance ze simulace
result_means = zeros(1,7); 
result_varis = zeros(1,7); 

for i = 1:length(chyby_odhadu)
    chyba_odhadu = chyby_odhadu{i}; 
    result_means(1, i) = mean(chyba_odhadu(1,:)); 
    result_varis(1, i) = var(chyba_odhadu(1,:)); 
end
