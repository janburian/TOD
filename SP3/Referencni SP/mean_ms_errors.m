function [output] = mean_ms_errors(e)

sim_count = length(e);
sim_length = length(e{1});

output = zeros(2, sim_length);

for i=1:sim_length
    e_mean = [0 0]';
    for j=1:sim_count
        e_mean = e_mean + e{j}(:, i);
    end
    e_mean = e_mean / sim_count;
    output(:, i) = e_mean;
end
end

