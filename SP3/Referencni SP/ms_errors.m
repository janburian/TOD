function [output] = ms_errors(vec1, vec2)
len = length(vec1(1, :));
if len ~= length(vec2(1, :))
    disp('ERR: Vectors must have same length.')
end
output = zeros(2, len);

for i=1:len
   output(:, i) = ( vec1(:, i) - vec2(:, i) ).^2; 
end

end

