
function [rme, b] = demarket(r, mkt, b)
    
    if ~exist('b', 'var') % optional 
        % compute market beta 
        rhs = [ones(size(mkt,1),1) mkt];
        b = rhs \ r; b = b(2, :);
    end
    
    % de-market
    rme = r - mkt*b;
end

