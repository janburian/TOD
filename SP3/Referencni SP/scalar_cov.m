function [output] = scalar_cov(vec1, vec2)

len = length(vec1);
if len ~= length(vec2)
    disp('ERR: Vectors must have same length.')
end

E_vec1 = sum(vec1)/len;
E_vec2 = sum(vec2)/len;

output = 0;

for i=1:length(vec1)
    output = output + (vec1(i) - E_vec1)*(vec2(i) - E_vec2);
end
output = output/len;

end

