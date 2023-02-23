clc
clear all
close all 

%% Zadany system 
%x_{k+1} = F*x_k + G*w_k, w_k ~ N(0,Q) 
%z_{k} = H*x_k + v_k, v_k ~ N(0,R)
     
%% Vygenerovani 1000 realizaci nahodneho vektoru z 
N = 1000; 
pocet_mereni = 4; 

Z = randn(N, pocet_mereni); 

%% Vypoctene estimatory 
est_z0_z1 = [1 0; 
             -1 1]; 
         
est_z0_z2 = [1 0; 
             -0.5 0.5];
         
est_z1_z2 = [2 -1;
             -1 1]; 
         
est_z0_z1_z2 = [76/91 30/91 -15/91; 
                -185/364 3/182 179/364]; 
         
est_z0_z1_z2_z3 = [3569/5025 652/1675 152/1675 -956/5025; 
                -24559/75375 -1772/25125 2978/25125 20941/75375]; 
            
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
figure
scatter(res_z0_z1(1,:), res_z0_z1(2,:), '.'); 

figure
scatter(res_z0_z2(1,:), res_z0_z2(2,:), '.'); 

figure
scatter(res_z1_z2(1,:), res_z1_z2(2,:), '.'); 

figure
scatter(res_z0_z1_z2(1,:), res_z0_z1_z2(2,:), '.'); 

figure
scatter(res_z0_z1_z2(1,:), res_z0_z1_z2(2,:), '.'); 

%% Odhad hustot jednotlivych odhadu polohy pomoci normalizovanych histogramu


