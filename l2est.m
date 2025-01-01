% L2 shrinkage (our main formula)
function [b, params, se] = l2est(X, y, params, compute_errors)
    if nargin < 4, compute_errors = false; end
    l = params.L2pen;

    if compute_errors 
        Xinv = inv(X + l*eye(size(X,1)));

        b = Xinv*y; 
        se = sqrt(1/params.T * diag(Xinv));
    else
        % solve a system of linear equation instead if errors are not needed
        b = (X + l*eye(size(X,1))) \ y; 

        se = nan(size(X, 1), 1);
    end
end

%         b = l*((X + 0*eye(size(X,1))) \ y);
%         b = (X + l*X^-1) \ y;
%         b = (X'*X + l*eye(size(X,1))) \ X'*y;
%         b = (X'*X + l*X) \ X'*y;
