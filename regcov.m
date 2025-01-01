function [X] = regcov(r)

    X = cov(r);

    % Covariance regularization (with flat Wishart prior)
    [T, n] = size(r);
    a = n / (n + T);
    X = a*trace(X)/n*eye(n) + (1-a)*X;
end