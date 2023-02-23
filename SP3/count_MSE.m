function [result] = count_MSE(cell_MSE)
%COUNT_MSE Summary of this function goes here
    matrix = [cell_MSE{:}]; % prevod na matici 
    MSE_poloha = matrix(:, 1:2:end); % poloha - liche sloupce
    MSE_rychlost = matrix(:, 2:2:end); % rychlost - sude sloupce
    
    dim = size(MSE_poloha);
    pocet_kroku = dim(1); 
    result = zeros(pocet_kroku, 2); 

    % Vysledny vypocet MSE
    for k = 1:pocet_kroku
        poloha = MSE_poloha(k, :); % hodnota polohy v case k pro vsechny trajektorie
        rychlost = MSE_rychlost(k, :); % hodnota rychlosti v case k pro vsechny trajektorie
        
        result(k, :) = [mean(poloha) mean(rychlost)]; 
    end

end

