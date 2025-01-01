% Please cite as: Kozak, Nagel, and Santosh (2019): ``Shrinking the Cross 
% Section'', Journal of Financial Economics.
%
% Copyright 2019, Serhiy Kozak, Stefan Nagel, and Shri Santosh. 
%
% First Version: 2017-04-20. Serhiy Kozak.
% This  Version: 2019-08-08. Serhiy Kozak.


%% main execution file; calls up scs_estimate.m for different test assets
clear; clear global; clc; close all;
tic

%% options 
% use daily returns throughout
daily = true;

% compute interactions
interactions = false;

% rotate returns into PC space 
rotate_PC = false;

% withhold test sample for OOS tests
withhold_test_sample = false;

% type of portfolios to use as primary data (pick any)
dataprovider = 'anom';
% dataprovider = 'ff25';

% sample dates
t0 = datenum('01-Jul-1963');
tN = datenum('31-Dec-2017');
oos_test_date = '01JAN2005';

% current run
run_folder = [upper(datestr(today, 'ddmmmyy')) '/'];

       
%% paths
projpath = './'
datapath = [projpath 'Data/']

instrpath = [datapath 'instruments/'];

%% initialize
% daily file suffix
if daily
    freq = 252;
    suffix = '_d';
    date_fmt = 'mm/dd/yyyy';
else
    freq = 12;
    suffix = '';
    date_fmt = 'mm/yyyy';
end

% Fix random number generation
rng default;

%% default estimation parameters
default_params = struct( ...
    'gridsize', 100, ... % size of the grid
    'contour_levelstep', 0.01, ... % L1L2 contour map plot: level step of contour lines 
    'objective', 'CSR2', ... % objective used in the paper
    'rotate_PC', false, ... % pre-rotate assets into the PC space
    'devol_unconditionally', false, ... % should be true unless constant-leveraged (already scaled) portfolios are used
    'kfold', 3, ... % 3-fold cross-validation
    'plot_dof', true, ... % show plots
    'plot_coefpaths', true, ...
    'plot_objective', true);
default_params.fig_options.fig_sizes = {'width=half'}; 
% default_params.fig_options.fig_sizes = {'width=full', 'width=half'}; % produce pictures in two sizes
default_params.fig_options.close_after_print = true; % close figure after saving

% load FF factors and save for future use (ALWAYS RUN THIS PIECE OF CODE)
[dd, re, ~] = load_ff_anomalies(datapath, daily, t0, tN);
global ff_factors ff_factors_dates; ff_factors = re; ff_factors_dates = dd;


% parameters
p = default_params;

if interactions % 2-fold CV for interactions
    p.kfold = 2;
else % fine grid for runs with no interactions
    p.gridsize = 100;
end

if withhold_test_sample % tests using withheld sample
    p.oos_test_date = oos_test_date;
end

%% process original ff25 portfolios if requested
if strcmp(dataprovider, 'ff25')
    if ~interactions % no interactions for FF25 -- simply ignore the request
        %%
        [dd, re, mkt, DATA, labels] = load_ff25(datapath, daily, 0, tN); % t0 -> 0
        p.smb = DATA.SMB; p.hml = DATA.HML; 
        p.L2_table_rows = 10; % 10 rows for FF25
        p.table_L2coefs_posneg_sort = ~rotate_PC; % equal number of positive & negative coefs to show (5)
        p.table_L2coefs_extra_space = rotate_PC; % add extra lspace for FF25 PCs (so that two tables have the same height)
        p.L2_sort_loc = 'OLS'; % sort labels at the rightmost point (~OLS solution)
        p.devol_unconditionally = true; % for FF25 we still need to de-vol since characteristics are not kept at const leverage
        p = SCS_L2est(dd, re, mkt, freq, labels, p);
    end
else
    %% Managed portfolios

    % find the appropriate file by mask
    fmask = [instrpath 'managed_portfolios_' dataprovider suffix '_*.csv'];
    flist = ls(fmask); 
    filename = [instrpath strip(flist(1,:))]; %assert(length(filename)>length(instrpath), ['File not found: ' fmask]);

    % parameters
    p.L1_truncPath = true; 

    if interactions % use interactions
        [dd, re, mkt, anomalies] = load_managed_portfolios(filename, daily, 0.2, {});
        p = SCS_L2est(dd, re, mkt, freq, anomalies, p);
    else % use only raw characteristics (no derived instruments)
        %% load data
        [dd, re, mkt, anomalies] = load_managed_portfolios(filename, daily, 0.2, {'rX_', 'r2_', 'r3_'});
        
        %% estimate
        p = SCS_L1L2est(dd, re, mkt, freq, anomalies, p);
    end
end

p.optimal_model_L2
