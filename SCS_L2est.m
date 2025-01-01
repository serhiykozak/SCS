% complete set of tests for a given portfolio set; this file gets called up by main.m

function [estimates] = SCS_L2est(dates, re, market, freq, anomalies, p) % p is a structure of parameters
%
% INPUTS
% - dates  : (T x 1) time series of numeric dates (in Matlab's format)
% - re     : (T x N) matrix of excess returns time series
% - market : (T x 1) matrix of market's excess returns time series 
% - freq   : number of obs per year
% - anomalies : cell array of anomaly names
% - p   : structure that contains extra arguments:
%     - 'gridsize'  : Default: 20. 
%                     Grid size for coefficient paths estimates.
%     - 'cvmethod'  : 'CV' (default), 'ssplit', 'bootstrp', 'AICc'. 
%                     Cross-validation method.
%     - 'kfold'     : Default: 4. 
%                     Number of folds in k-fold cross validation
%     - 'objective' : 'C-S R^2' (default), 'SSE', 'HJ dist'. 
%                     Objective function for cross validation.
%     - 'plot_dof'      : Default: false. Plot degrees of freedom.
%     - 'plot_coefpaths': Default: false. Plot L2 coefficient paths.
%     - 'plot_objective': Default: false. Plot model selection diagnostics.
%
%     - see the code below for more parameters and options...
%
% OUTPUTS
% - estimates : structure of estimates with following fields
%     - 'optimal_model': structure that contains the optimally selected model:
%          - 'coefficients': (N x 1) vector of coefficients.
%          - 'se'          : (N x 1) vector of coefficients' standard errors.
%          - 'objective'   : value of objective (OOS CS R^2 by default)
%          - 'dof'         : effective degrees of freedom (strength of
%                            regularization)
%          - 'kappa'       : optimalkappa (Root Expected SR^2)
%
%     - 'coeffPaths': (N x gridsize) matrix of coefficients.
%     - 'df'        : (1 x gridsize) vector of effective degrees of freedom
%                     at each grid point.
%     - 'objL2_IS'  : (gridsize x 1) vector of in-sample objective values
%                     at each grid point.
%     - 'objL2_OOS' : (gridsize x 1) vector of out-of-sample objective values
%                     at each grid point.
%
%
% DESCRIPTION: Computes the L2 shrinkage estimator of the SDF parameters
% based on the method in Kozak, Nagel, and Santosh (2019). 
% See the paper for further details.
%
% Please cite as: Kozak, Nagel, and Santosh (2019): ``Shrinking the Cross 
% Section'' Working Paper, University of Michigan.
%
% Copyright 2019, Serhiy Kozak, Stefan Nagel, and Shri Santosh. 
%
% First Version: 2017-04-20. Serhiy Kozak.
% This  Version: 2019-08-08. Serhiy Kozak.

%% validate arguments
assert(sum(sum(isnan(re),2)>0) == 0, 'Missing obs. found.')

%% assign default values

% general options
P = struct;
P.gridsize = 20; % number of points on the L2 grid in kappas
P.method = 'CV'; % the only method properly supported and backtested 
P.objective = 'CSR2'; % cross-validation objective
P.ignore_scale = false; % set this to true to allow rescaling fitted SDF in the OOS period
P.kfold = 5; % number of folds in K-fold cross-validation
P.oos_test_date = datestr(dates(end)); %'01JAN1995'; % set this to withhold part of the sample for OOS tests
P.freq = freq; % number of obs per year
P.rotate_PC = false; % rotate returns into PC space
P.demarket_conditionally = false; % de-market conditionally?
P.demarket_unconditionally = true; % de-market unconditionally, if demarket_conditionally=false?
P.devol_conditionally = false; % de-vol unconditionally if false
P.devol_unconditionally = true; % skip de-vol unconditionally if false

% figure options
P.plot_dof = false; % Plot degrees of freedom
P.plot_coefpaths = false; % Plot L2 coefficient paths
P.plot_objective = false; % Plot model selection diagnostics
P.line_width = 1.5; % figure options
P.font_size = 10; % default font size for legends
P.L2_max_legends = 20; % number of legends to plot on L2 coeffs paths plots
P.L2_sort_loc = 'opt'; % sort labels either at optimal kappa (opt) or at th epoint of no shrinkage ('OLS')
P.L1_log_scale = true; % % use log scale for (# of variables)
P.L2_log_scale = true; % use log scale for kappa
P.legend_loc = 'bestoutside';

% parse config and assign default values; options can be a string
p = parse_config(p, P);

% we usually maximize an objective (e.g., R^2), except for HJ-distance (GLS) and SSE
if any(strcmp(p.objective, {'GLS', 'SSE'}))
    optfunc = @min;
else
    optfunc = @max;
end

% user-friendly names for objectives to use in plots
mapObj = containers.Map({'CSR2', 'GLSR2', 'GLS', 'SRexpl', 'SSE', 'SR', 'MVU'}, ...
    {'Cross-sectional $R^2$', 'Cross-sectional GLS $R^2$', 'Residual $SR^2$', 'Explained SR', ...
    'SDF RMSE', 'Sharpe Ratio', 'Mean-variance utility'});
p.sObjective = mapObj(p.objective);

%% initialize; compute means, cov, SVD decomposition
% testing sample start date
tT0 = datenum(p.oos_test_date);

% train set is used for fitting AND cross-validation
idx_train = find(dates <= tT0);
% test set is optional and is used to full OOS testing of already fitted
% and cross-validated model. It is often empty (not used)
idx_test = find(dates > tT0);


mkt0 = market;

% de-market all excess returns 
if p.demarket_conditionally % conditionally
    demarket_ma_window = 3*freq; % use past 3 years to estimate betas
    r0 = demarketcond(re(idx_train,:), market(idx_train), demarket_ma_window);
    idx_train = idx_train(1+demarket_ma_window:end); % drop nans
elseif p.demarket_unconditionally % unconditionally
    [r_train, b_train] = demarket(re(idx_train,:), market(idx_train,:));
    r_test = demarket(re(idx_test,:), market(idx_test,:), b_train); % use betas estimated in the training sample
    r0 = [r_train; r_test];
else
    r0 = re;    
end
 
% de-vol all excess returns conditionally if requested
if p.devol_conditionally
    devol_ma_window = 22; % use past 22 days to estimate volatilities
    [r0, mkt0] = devolcond(r0, market, devol_ma_window);
    idx_train = idx_train(1+devol_ma_window:end); % drop nans
elseif p.devol_unconditionally % de-vol unconditionally
    % normalize so that all returns have standard deviation of the VW market
    r0 = r0./ repmat(nanstd(r0),size(dates,1),1) * nanstd(market);
end

% construct dates, mkt, and returns for train and test sets
dd = dates(idx_train);
dd_test = dates(idx_test);
mkt = mkt0(idx_train);
mkt_test = mkt0(idx_test);

r_train = r0(idx_train,:);
r_test = r0(idx_test,:);

% length of the training sample
[T, n] = size(r_train); 
p.T = T; p.n = n;


% rotate into PC space if requested and change file suffix
if p.rotate_PC
    [Q, ~] = svd(regcov(r_train)); % use training sample to form eigenvectors
    r_train = r_train*Q;
    r_test = r_test*Q;
    anomalies = strcat('PC',strtrim(cellstr(num2str((1:n)'))))';
end

% compute first and second moments
X = regcov(r_train); y = mean(r_train)';
X_test = regcov(r_test); y_test = mean(r_test)';

% maximum in-sample SR
w = y'*pinv(X); maxSR2 = freq*w*y;
% w_test = y_test'*pinv(X_test); maxSR2_test = freq*w_test*y_test

% precompute E-V decomposition
[Q, D] = svd(X);
X2 = Q * D.^.5 * Q';
d = diag(D);

% pre-compute pseudo inverses
tol = max(size(X)) * eps(norm(d,inf));
r1 = sum(d > tol)+1; Q1 = Q; s = d;
Q1(:,r1:end) = [];
s(r1:end) = [];
s2 = 1./sqrt(s(:)); s = 1./s(:); 
Xinv = (Q1.*s.')*Q1';
X2inv = (Q1.*s2.')*Q1';

% options
p.xlbl = 'Root Expected SR$^2$ (prior), $\kappa$';
p.Q = Q; 
p.d = d;
p.Xinv = Xinv;

%% Create L2 grid with flexible bounds
% functions to map L2pen <-> kappa
kappa2pen = @(kappa, T, X, p) p.freq*trace(X)/T ./ (kappa.^2);
% note that the sample size is shortened due to data with holding required for cross-validation

% find left and right limits
lr = 1:21; lm = 1;
z = nan(n, length(lr));
for i = lr 
%     params = p; params.L2pen = freq*trace(X)/(T/2) / 10^(i-lm);
    params = p; params.L2pen = kappa2pen(2^(i-lm), T, X, p); % define grid in kappas
    z(:, i) = l2est(X, y, params);
end
% coefficient stabilize when they do not change by more than 1% on average
% for penalties that differ by a factor of 10
x_rlim = lr(mean( abs((z(:,2:end)-z(:,1:end-1)))./(1+abs(z(:,1:end-1))) ) > 0.01)-lm;

% use the left and right points to define the support and create a finer
% grid on this support
% set minimum kappa at 0.01; estimate max above
x = logspace(log10(2^x_rlim(end)), log10(0.01), p.gridsize); % x variable for all plots: Root Expected SR^2 (kappa)
l = kappa2pen(x, T, X, p); % L2 penalty corresponding to kappa for estimates based of FULL sample
lCV = l / (1-1/p.kfold); % adjust L2 penalty due to shorter sample used in CV estionation (same as above)
nl = length(l);

%% Estimate the L2 model
params = p; % for GLS models cross_validate() can pre-cache matrix inverses 
            % between the calls -- avoid parfor
% estimate the L2 problem
phi = nan(n, nl); se = phi;
objL2 = nan(nl, 4); 
MVE = cell(nl, 1); % list of MVE portfolio returns for each level of regularization
for i = 1:nl 
    % estimate parameters at each grid point
    params.L2pen = l(i); % penalty for estimates based on full sample
    [phi(:, i), ~, se(:, i)] = l2est(X, y, params, true);

    % cross validate estimated parameters
    params.L2pen = lCV(i); % penalty for estimates based on shorter CV sample
    [objL2(i, :), params, objL2_folds_] = cross_validate(@l2est, dd, r_train, params);
    objL2_folds(i, :) = objL2_folds_(:, 2); % save OOS R2 estimates for all folds    

    % store OOS MVE portfolios for each CV run
    MVE{i} = params.cv_MVE;
end
cv_idx_test = params.cv_idx_test; % CV test sample indices

% effective degrees of freedom
df = sum(d./(d + l)); % based on full sample penalty (same as eigenvalues)

% optimal L2 model
[objL2opt, iL2opt] = optfunc(objL2(:, 2));
bL2 = phi(:, iL2opt);
p.bL2 = bL2;
p.R2oos = objL2opt; 
L2optKappa = x(iL2opt);

% MVE portfolios for each fold at the optimal level of shrinkage [flatten into single time-series]
MVEopt = MVE{iL2opt};

% return coefficients paths, degrees of freedom, and objective's value
p.coeffsPaths = phi;
% p.df = df; 
p.objL2_IS = objL2(:, 1);
p.objL2_OOS = objL2(:, 2);
p.optimal_model_L2.coefficients = bL2;
% p.optimal_model.se = se(:, iL2opt);
p.optimal_model_L2.objective = objL2opt;
% p.optimal_model.dof = df(iL2opt);
p.optimal_model_L2.kappa = L2optKappa;
z = cell2mat(MVEopt'); % [MVEopt{1}; MVEopt{2}; MVEopt{3}; MVEopt{4}; MVEopt{5}];
p.optimal_model_L2.SR = mean(z) / std(z) * sqrt(p.freq);

estimates = p;


%% df <-> kappa plot
if p.plot_dof % plot degrees of freedom
    plot_dof(df, x, p);
end

%% SDF 2nd moment constraint (L2) coefficients 
if p.plot_coefpaths 
    % plot coefficients
    plot_L2coefpaths(x, phi, iL2opt, anomalies, 'SDF Coefficient, $b$', p);
    % plot t-stats
    plot_L2coefpaths(x, phi./se, iL2opt, anomalies, 'SDF Coefficient $t$-statistic', p);
end
    
%% L2 Cross-Validation/BIC plot
if p.plot_objective
    plot_L2cv(x, objL2, p);
end

%% output table with coefficient & tstats estimates
table_L2coefs(phi(:, iL2opt), se(:, iL2opt), anomalies, p)


end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% degrees of freedom <-> kappa plot
function [] = plot_dof(df, x, p)
%%
    % open new figure
    figure(); 

    % plot
    plot(x, df, 'LineWidth', p.line_width)
    if p.L1_log_scale
        set(gca,'yscale','log'); set(gca,'YTick',get(gca,'YTick')+1e-12) % adding a small constant changes ticks to integer numbers rather than 10^x powers
    end
    if p.L2_log_scale
        set(gca,'xscale','log'); set(gca,'XTick',get(gca,'XTick')+1e-12) % adding a small constant changes ticks to numbers instead of 10^k.
    end
    
    xlabel(p.xlbl, 'Interpreter', 'latex')
    ylabel('Effective degrees of freedom', 'Interpreter', 'latex')
    grid on

    set(gca, 'XLim', [min(x) max(x)]);
end

%L2 coefficients paths plot
function [] = plot_L2coefpaths(x, phi, iL2opt, anomalies, ylbl, p)
%%
    switch p.L2_sort_loc
        case 'opt'
            iSortLoc = iL2opt;
        case 'OLS'
            iSortLoc = 1;
        otherwise 
            error('Unknown option');
    end
    % if more variables than L2_max_legends, show only the largest (abs. value)
    % otherwise show all and order from top to bottom
    if p.n > p.L2_max_legends
        [~, I] = sort(abs(phi(:, iSortLoc)), 'descend');
    else
        [~, I] = sort(phi(:, iSortLoc), 'descend');
    end

    % open new figure
    figure()
    hold on

    % plot
    plot(x, phi(I,:)', 'LineWidth', p.line_width); 
    if p.L2_log_scale
        set(gca,'xscale','log'); set(gca,'XTick',get(gca,'XTick')+1e-16)
    end
    grid on; 
    xlabel(p.xlbl, 'Interpreter', 'latex'); 
    ylabel(ylbl, 'Interpreter', 'latex')

    % add top L2_max_legends legends
    idx = I(1:min(p.L2_max_legends,length(I)));
    legend(anomnames(anomalies(idx)), 'Location', p.legend_loc, ...
        'FontSize', p.font_size, 'Interpreter', 'latex'); legend('boxoff')
    
    % plot a dashed line at the optimal level of regularization
    plot(x(iL2opt)*[1 1], [min(min(phi)) max(max(phi))], '--k')

    set(gca, 'XLim', [min(x) max(x)]);
end

% plot SSE/objective & BIC as a function of degrees of freedom
function [] = plot_L2cv(x, objL2, p)
%%
    % plot
    figure(); 
    plot(x, objL2(:, 1), '--', 'LineWidth', p.line_width); hold on % IS
    plot(x, objL2(:, 2), '-', 'LineWidth', p.line_width); % OOS
    
    if p.L2_log_scale
        set(gca,'xscale','log'); set(gca,'XTick',get(gca,'XTick')+1e-16)
    end
    xlabel(p.xlbl, 'Interpreter', 'latex')
    ylabel(['IS/OOS ' p.sObjective], 'Interpreter', 'latex');

    legends = {'In-sample', ['OOS ' p.method]};
    
    % plot se
    co = get(gca, 'ColorOrder');
    plot(x, objL2(:, 2)+objL2(:, 4), ':', 'Color', co(2,:), 'LineWidth', 1); 
    plot(x, objL2(:, 2)-objL2(:, 4), ':', 'Color', co(2,:), 'LineWidth', 1); 

    legend([legends ['OOS ' p.method ' \pm 1 s.e.']], ...
        'Location', 'NorthWest'); legend('boxoff'); grid on; 

    grid on
    ylim([0 max(0.1,min(10,2*max(objL2(:, 2))))])
    xlim([min(x) 2])
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% a table listing the largest coefficients and t-stats
function [] = table_L2coefs(phi, se, anomalies, p)
    nrows = 10; % number of rows in the table to show
    
    % t-stats
    tstats = phi./se;
    
    % by absolute tstats
    [~, I] = sort(abs(tstats), 'descend');

    % show only L2_table_rows items
    idx = I(1:nrows);   
   
    % output table
    tbl = table(anomdescr(anomalies(idx)), phi(idx), abs(tstats(idx)), 'VariableNames', {'Portfolio', 'b', 't_stat'});
    disp(tbl)
end
