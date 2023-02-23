function [stopy_vektor] = count_traces(MSE_teo_cell)
%COUNT_TRACES Summary of this function goes here
%   Detailed explanation goes here
stopy_vektor = zeros(26,1); 
    for j = 1:length(MSE_teo_cell)
        matrix = MSE_teo_cell{j}; 
        stopa = trace(matrix); 
        stopy_vektor(j, 1) = stopa; 
    end
end

