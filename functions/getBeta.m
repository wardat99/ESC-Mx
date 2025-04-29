function [beta] = getBeta(n3)


% Initialize beta with the first three fixed values
beta = [0.5; 5; 100];

% If n3 > 3, add more values increasing by 10
if n3 > 3
    extra_values = (100 + 10) + (0:(n3-4)) * 10;  % Start from 110, step 10
    beta = [beta; extra_values'];
end


end