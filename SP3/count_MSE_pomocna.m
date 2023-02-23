function [error] = count_MSE_pomocna(trajektorie_vek, odhad_vek)
%COUNT_MSE Summary of this function goes here
    dimenze_vek = size(trajektorie_vek); 
    error = zeros(dimenze_vek(1), dimenze_vek(2)); 
    
    for j = 1:length(trajektorie_vek)
            error(j, :) = (trajektorie_vek(j, :) - odhad_vek(j,:)).^2;  
    end 
end

