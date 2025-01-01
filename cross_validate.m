function [obj, params, obj_folds] = cross_validate(FUN, dates, r, params)
% INPUTS
% - FUN    : Handle to a function which estimates model paramaters
% - date   : (T x 1) vector of dates
% - r      : (T x N) matrix of returns 
% - params : structure that contains extra arguments:
%     - 'method'    : 'CV' (default), 'ssplit', 'bootstrap'. 
%                     Cross-validation method.
%     - 'objective' : 'SSE', 'GLS', 'CSR2', 'GLSR2', 'SR', 'MVU'. 
%                     Objective function for cross validation.
%     - 'ignore_scale' : allows for an additional degree of freedom when
%                        calculating IS/OOS goodness-of-fit statistics by
%                        re-scaling the MVE portfolio for best fit (regress
%                        mean returns on MVE portfolio return and scale up
%                        MVE portfolio using the estimated coefficient).
%                        Default: false.
%     - 'kfold'     : Number of folds in k-fold cross validation
%     - 'splitdate' : Date of sample split in sample split validation
%
% OUTPUTS
% - obj: (1 x 2) IS and OOS values of the estimated objective function
% - params: returns back the params structure. If cross_validate() is
%           called multiple times in a row, it may use params to cache
%           some calculations in between those calls (currently matrix
%           inverses are cached for GLS objectives)
%
% DESCRIPTION: Compute IS/OOS values of the objective function based on the
%              FUN function. Implements multiple objectives and validation
%              methods.
%
% Please cite as: Kozak, Nagel, and Santosh (2019): ``Shrinking the Cross 
% Section'', Journal of Financial Economics.
%
% Copyright 2019, Serhiy Kozak, Stefan Nagel, and Shri Santosh. 
%
% First Version: 2017-07-13. Serhiy Kozak.
% This  Version: 2017-07-13. Serhiy Kozak.


    narginchk(3, 5);
    if ~isa(FUN,'function_handle')
        error(message('MATLAB:cross_validate:funArgNotHandle'));
    end
    
    % select requested method
    if ~isfield(params, 'method')
        cross_validate_def_handler = @cross_validate_cv_handler;
    else
        map_cv_method = containers.Map({'CV', 'ssplit', 'bootstrap'}, ...
            {@cross_validate_cv_handler, @cross_validate_ssplit_handler, @cross_validate_bootstrap_handler});

        cross_validate_def_handler = map_cv_method(params.method);
    end
    
    % execute selected method
    params.dd = dates;
    params.ret = r;
    params.fun = FUN;
    [obj, params, obj_folds] = cross_validate_def_handler(params);
end

% Sample split method
function [obj, params] = cross_validate_ssplit_handler(params)
    if isfield(params, 'splitdate')
        sd = params.splitdate;
    else
        sd = '01JAN2000';
    end

    tT0 = datenum(sd);
    idx_test = find(params.dd>=tT0);

    [obj, params] = bootstrp_handler(idx_test, params);
end

% Bootstrap method
function [obj, params] = cross_validate_bootstrap_handler(params)
    error('Not implemented properly!')
    
    idx = (1:length(params.r))';
    opt = statset('UseParallel', false);
    bootstat = bootstrp(params.niter, @bootstrp_handler, idx, params, 'Options', opt); % CHECK THIS
    
    obj = mean(bootstat, 1);
end

% Cross-validation method   
function [obj, params, obj_folds] = cross_validate_cv_handler(params)
    if isfield(params, 'kfold')
        k = params.kfold;
    else
        k = 2;
    end
    
    % We use contiguous partitioning for financial data tue to their high
    % persistence. Using random sampling partitioning produces highly
    % correlated samples
    cv = cvpartition_contiguous(size(params.ret,1), k); % contiguous partition
%     cv = cvpartition(size(params.ret,1), 'KFold', k); % random partition
    
    % compute IS/OOS statistics for each partition
    obj = nan(k, 2);
    for i = 1:k
        idx_test = cv{i}; % contiguous partition
        params.cv_idx_test{i} = idx_test; % store test indices
%         idx = find(training(cv, i)); % random partition
        
        params.cv_iteration = i;
        
        [obj(i, :), params] = bootstrp_handler(idx_test, params);
    end
    
    % return estimates for each fold
    obj_folds = obj;
    
    % average IS/OOS statistics across folds
    obj = [mean(obj, 1) std(obj, 1)/sqrt(k)];
    
%     % take square root for SRexpl
%     if strcmp(params.objective, 'SRexpl')
%         obj = sqrt(max(0,obj));
%     end
end


% This function executes the logic of of the problem
function [res, params] = bootstrp_handler(idx_test, params)

    % pick requested objective
    if isfield(params, 'objective')
        map_bootstrp_obj = containers.Map({'SSE', 'GLS', 'CSR2', 'GLSR2', 'SRexpl', 'SR', 'MVU'}, ...
            {@bootstrp_obj_SSE, @bootstrp_obj_HJdist, @bootstrp_obj_CSR2, ...
            @bootstrp_obj_GLSR2, @bootstrp_obj_SRexpl, @bootstrp_obj_SR, @bootstrp_obj_MVutil});
        
        def_bootstrp_obj = map_bootstrp_obj(params.objective);
    else
        def_bootstrp_obj = @bootstrp_obj_SSE;
    end

    % initialize parameters
    ret = params.ret;
    FUN = params.fun;
    
    % construct test sample indices
    n = size(ret,1);
    idx = (1:n)'; % training indices
    idx(idx_test) = [];
    n_test = size(idx_test,1);

    invX = nan; invX_test = nan;
    res = [nan nan];
    if n_test > 0 % test sample non-empty
        % IS/OOS returns
        r = ret(idx, :);
        r_test = ret(idx_test, :);

        % if no cache is available, calculate means and cov and cache them
        % for faster execution in case cross_validate() is called multiple times
        if ~isfield(params, 'cv_cache') || length(params.cv_cache) < params.cv_iteration
            % compute means and covariances
            cvdata = struct();
            cvdata.X = regcov(r); cvdata.y = mean(r)';
            cvdata.X_test = regcov(r_test); cvdata.y_test = mean(r_test)';
            if any(strcmp(params.objective, {'GLS', 'GLSR2', 'SRexpl'}))
                cvdata.invX = pinv(cvdata.X);
                cvdata.invX_test = pinv(cvdata.X_test);
            end
            
            % store in cache
            params.cv_cache{params.cv_iteration} = cvdata;
        end
        
        % read from cache 
        cvdata = params.cv_cache{params.cv_iteration};
        X = cvdata.X; y = cvdata.y;
        X_test = cvdata.X_test; y_test = cvdata.y_test;
        if any(strcmp(params.objective, {'GLS', 'GLSR2', 'SRexpl'}))
            invX = cvdata.invX;
            invX_test = cvdata.invX_test;
        end
        
        % estimate a model
        [phi, params] = FUN(X, y, params); 
        
        if ~isfield(params, 'cache_run') || ~params.cache_run

            % store coefficient estimates for each CV run
            params.cv_phi{params.cv_iteration} = phi;
            % store OOS MVE portfolios estimates for each run 
            params.cv_MVE{params.cv_iteration} = r_test*phi;

            % estimate the IS/OOS MVE portfolio
            fact = X*phi;
            fact_test = X_test*phi;

            % undo scaing if requested
            if params.ignore_scale
                % fit OLS
                b = fact \ y;
                b_test = fact_test \ y_test;  
            else
                % no rescaling
                b = 1;
                b_test = 1;
            end

            % return IS & OOS summary statistics
            res = [ def_bootstrp_obj(fact*b, y, invX, phi, r, params) ... % IS 
                    def_bootstrp_obj(fact_test*b_test, y_test, invX_test, phi, r_test, params) ... % OOS
                  ];
        end
    end    
end

% GLS/HJdist objective
function [obj] = bootstrp_obj_HJdist(y_hat, y, invX, phi, r, params)
    obj = sqrt((y_hat - y)'*invX*(y_hat - y) * params.freq); % unexplained squared Sharpe Ratio
end

% SSE objective
function [obj] = bootstrp_obj_SSE(y_hat, y, invX, phi, r, params)
    obj = params.freq * sqrt((y_hat - y)'*(y_hat - y) / size(y,1));
end

% Cross-sectional R^2 objective under the null (model) that all alphas are zero
function [obj] = bootstrp_obj_CSR2(y_hat, y, invX, phi, r, params)
    obj = 1 - (y_hat - y)'*(y_hat - y) / (y'*y); 
end

% Predictive R^2 objective (panel)
function [obj] = bootstrp_obj_predR2(y_hat, y, invX, phi, r, params)
    e = r - y_hat'; e = e(:);
    obj = 1 - e'*e / (r(:)'*r(:)); 
end

% Cross-sectional GLS R^2 objective under the null (model) that all alphas are zero
function [obj] = bootstrp_obj_GLSR2(y_hat, y, invX, phi, r, params)
    obj = 1 - (y_hat - y)'*invX*(y_hat - y) / (y'*invX*y);
end

% Sharpe ratio explained (diff in HJ dist of our and risk-neutral model)
% for small PCs, y_hat_PC ~ 0, so small PCs cancel out...
% see Barillas, Shanken 2017 p.7
function [obj] = bootstrp_obj_SRexpl(y_hat, y, invX, phi, r, params)
    obj = ((y'*invX*y) - (y_hat - y)'*invX*(y_hat - y)) * params.freq; 
end

% Sharpe Ratio objective
function [obj] = bootstrp_obj_SR(y_hat, y, invX, phi, r, params)
%     sd = tsmovavg((r.^2)', 's', params.freq/2)';
%     r = (r*phi)./sd * mean(std(r));
    obj = sqrt(params.freq)*nanmean(r*phi)/nanstd(r*phi);
end

% Mean-Variance Utility objective
function [obj] = bootstrp_obj_MVutil(y_hat, y, invX, phi, r, params)
    pret = 1+r*phi; 
    if any(pret < 0)
        obj = -1e9;
    else
        assert(all(pret > 0));
        pret = log(pret);
        obj = params.freq * (mean(pret) - 2*0.5*var(pret)); % 2* is the Jensen's correction (log vs levels)
    end
%     obj = params.freq * (mean(pret/params.gamma) - params.gamma/2*var(pret/params.gamma));
end
