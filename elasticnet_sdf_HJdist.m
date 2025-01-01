%  Front end for LARS-EN algorithm in SPASM package to implement double SDF
%  shrinkage (L1 and L2) in Kozak, Nagel, Santosh (2017). 
%  GLS objective that minimizes HJ distance between a target and
%  unrestricted SDFs: min (Y - X*b)' * S^-1 * (Y - X*b) s.t. b'b<=B and L1 penalty
%  Importantly, none of the variables are centered or standardized. Pricing equation in KNS (2017)
%  should hold without intercept and using original (not rescaled) predictors.
%  
%
% Copyright 2017, Serhiy Kozak. 
%
% First Version: 2017-07-21. Serhiy Kozak.
% This  Version: 2017-07-21. Serhiy Kozak.
%


function [b, params] = elasticnet_sdf_HJdist(X, y, params) % X2, X2inv, delta, stop, storepath, verbose)
%ELASTICNET Regularization and variable selection for regression via the
%   Elastic Net.
%
%   BETA = ELASTICNET_SDF(X, Y, DELTA) performs Elastic Net [1] regression on
%   the variables in X to approximate the response y. The regularization parameter DELTA
%   specifies the weight on the L2 penalty on the regression coefficients.
%   A positive value of DELTA makes analyses of p>n datasets possible and
%   encourages grouping of correlated variables. Setting DELTA to zero
%   yields the LASSO solution. BETA contains the Elastic Net solutions for
%   all values of the L1 regularization parameter, starting from a large
%   value corresponding to all beta coefficients equal to zero, to 0,
%   corresponding to the ordinary least squares solution.
%   BETA = ELASTICNET_SDF(X, Y, DELTA, STOP) with nonzero STOP yields Elastic
%   Net solutions with early stopping at a particular value of the L1
%   regularization parameter. If STOP is negative, STOP is an integer that
%   determines the desired number of non-zero variables. If STOP is
%   positive, it corresponds to an upper bound on the L1-norm of the BETA
%   coefficients. Setting STOP to zero (default) yields the entire
%   regularization path.
%   BETA = ELASTICNET_SDF(X, Y, DELTA, STOP, STOREPATH) with STOREPATH set to
%   true (default) will return the entire regularization path from zero
%   active variables (a high value of delta) to the point where the STOP
%   criterion is met, or the least squares solution in reached when STOP =
%   0. Setting STOREPATH to false will yield the single solution at which
%   the STOP criterion is met, thus saving precious computer resources.
%   BETA = ELASTICNET_SDF(X, Y, STOP, STOREPATH, VERBOSE) with VERBOSE set to
%   true (default is false) will print the adding and subtracting of
%   variables as LASSO solutions are found.
%   [BETA INFO] = ELASTICNET_SDF(...) returns a structure INFO containing
%   various useful measures. If the entire solution path has been
%   calculated, INFO containts the goodness-of-fit estimates AIC (Akaike's
%   Information Criterion), BIC (Bayesian Information Criterion), Mallow's
%   Cp statistic, and the number of degrees of freedom at each step. It
%   also includes the L1 penalty constraints s and lambda of which the
%   former represents the size of the L1 constraint defined in the range
%   [0,1], and the latter is lambda in the forumlation beta = arg min ||y
%   - X*beta||^2 + delta*||beta||^2 + lambda*||beta||_1. If stop <> 0
%   and a single Elastic Net solution is returned, INFO will only contain
%   the lambda and s value at the solution along with the number of
%   degrees of freedom. INFO also includes the number of steps made to
%   compute the solution, including the first step where beta is the zero
%   vector.
%
%   The algorithm is a variant of the LARS-EN algorithm [1].
%
%   References
%   -------
%   [1] H. Zou and T. Hastie. Regularization and variable selection via the
%   elastic net. J. Royal Stat. Soc. B. 67(2):301-320, 2005.
%

%% Input checking
narginchk(3, 3);


%% check whether the algorithm needs to be executed or parameters are already stored
if isfield(params, 'use_precomputed') 
    % do not execute the algorithm; simply read off stored coefficients
    bpath = params.bpath{params.cv_iteration}; % fint current CV iteration
    idx = find(sum(bpath~=0) >= abs(params.stop),1); % find coefficients
    b = bpath(:, idx);
    return
end


%% defaults
P = struct();
P.cv_iteration = 1; % current iteration (only for cross-validation)
P.L2pen = 0;
P.stop = 0;
P.storepath = false;
P.verbose = false;

% parse config and assign default values; options can be a string
p = parse_config(params, P);

assert(p.L2pen >= 0, 'L2 penalty must be non-negative.');

% parameters
[~, k] = size(X);
delta = p.L2pen;
b = [];

%% pre-compute matrix inverses and powers (these can be potentially precomputed by L2 method)
if isfield(p, 'elasticnet_cache') && length(p.elasticnet_cache)>=p.cv_iteration
    % read from cache
    cache = p.elasticnet_cache{p.cv_iteration};

    X1 = cache.X1;
    y1 = cache.y1;
else
    [Q, D] = svd(X);
    X2 = Q * D.^.5 * Q';
    d = diag(D);

    % pseudo inverses
    tol = max(size(X)) * eps(norm(d,inf));
    r1 = sum(d > tol)+1; Q1 = Q; s = d;
    Q1(:,r1:end) = []; s(r1:end) = [];
    s2 = 1./sqrt(s(:)); 
    X2inv = (Q1.*s2.')*Q1';

    % compute variables
    X1 = X2; % modified X
    y1 = X2inv*y; % modified y
    
    % save to cache
    cache = struct;
    cache.X1 = X1; 
    cache.y1 = y1; 
    params.elasticnet_cache{p.cv_iteration} = cache;
end

%% quit if pre-cache run is requested
if isfield(params, 'cache_run') && params.cache_run
    return
end

%% Calculate Elastic Net solutions with the LARS-EN algorithm

% temporarily treat warning as errors (to catch them in try..catch)
warnState(1) = warning('error', 'MATLAB:singularMatrix');
warnState(2) = warning('error', 'MATLAB:illConditionedMatrix');
warnState(3) = warning('error', 'MATLAB:nearlySingularMatrix');

% execute larsen() algorithm on modified variables (X1, y1)
bpath = larsen(X1, y1, delta, p.stop, [], p.storepath, p.verbose);
% bpath = larsenc_mex(X1, y1, delta, p.stop, Gram, p.storepath, p.verbose);

% adjust to avoid double shrinkage (non-naïve Elastic Net solution)
% bpath = (bpath~=0).*repmat(bpath(:, end), 1, size(bpath,2)); % undo L1 shrinkage
% bpath = ( 1 + delta/mean(d) ) * bpath;
% bpath = (1 + delta/mean(diag(X))) * bpath;
% % intelligently undo level shrinkage
% L1LevelShrinkage = 1-sum(abs(bpath))/sum(abs(bpath(:,end)));
% bpath = repmat(1 + delta/mean(d)*L1LevelShrinkage,size(bpath,1),1) .* bpath;

% save full coefs paths under current CV iteration (to reuse for future calls)
params.bpath{p.cv_iteration} = bpath;

% get coefficient values for a given L1 penalty (# of non-zero variables)
b = bpath(:, end);

% restore warnings
warning(warnState);

end

%%
%     % R implementation
%     tempfolder = ['m:/temp/R/' num2str(randi(1e9)) '/'];
%     mkdir(tempfolder);
%     csvwrite([tempfolder 'X.csv'],X1);
%     csvwrite([tempfolder 'y.csv'],y1);
% 
%     system(['"C:\Program Files\R\R-3.3.2\bin\Rscript.exe" elasticnet.r ' tempfolder ' ' num2str(delta)]);
%     bpath = csvread([tempfolder 'b.csv'])';
%     rmdir(tempfolder, 's')

% options = struct('alpha', 1, 'intr', 0, 'standardize', 0);
% options.lambda = [delta];
% tic; fit = glmnet(X1, y1, 'gaussian', options); toc; fit


% bpath = larsenc_mex(X1, y1, delta, stop, Gram, storepath, verbose);
