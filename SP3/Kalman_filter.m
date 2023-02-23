function [Kalman_result, inovace_KF] = Kalman_filter(trajektorie_cell, mereni_cell, pocet_kroku, F, G, H, Q, R, x0_cov, x0_mean)
%KALMAN_FILTER Summary of this function goes here
%   Detailed explanation goes here
    Kalman_result = cell(1,1000);
    vektor_stavu_Kalman = zeros(pocet_kroku+1, 2); 
    inovace_KF = zeros(pocet_kroku+1, 1000); 

    for i = 1:length(trajektorie_cell)
        predikce_P = x0_cov;
        predikce_x = x0_mean;

        mereni_vektor = mereni_cell{i}; 
        trajektorie = trajektorie_cell{i}; 

        mereni = mereni_vektor(1); 

        inovace_KF(1, i) =  mereni - predikce_x(1);
        korekce_P = predikce_P;
        korekce_x = predikce_x;

        vektor_stavu_Kalman(1,1:2) = korekce_x';    

        for k = 2:length(trajektorie)
            predikce_P = F * korekce_P * F' + G * Q * G';
            predikce_x = F * korekce_x; 

            mereni = mereni_vektor(k); 
            inovace_KF(k, i) = mereni - predikce_x(1);
            K = predikce_P * H' * (H * predikce_P * H' + R)^(-1);
            korekce_P = (eye(2) - K * H) * predikce_P * (eye(2) - K * H)' + K * R * K';
            korekce_x = predikce_x + K * (mereni - H * predikce_x);

            vektor_stavu_Kalman(k, 1:2) = korekce_x;
        end
        Kalman_result{i} = vektor_stavu_Kalman; 
    end
end

