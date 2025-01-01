function [indices] = cvpartition_contiguous(n, k)
    s = floor(n/k);
    indices = cell(k, 1);
    for i = 1:k-1
        indices{i} = s*(i-1)+1:s*i;
    end
    indices{k} = s*(k-1)+1:n;
end

