% Purpose: load a list of managed portfolios from SAS ouput
% 
% Input parameters:
%   datapath - path to the CSV file
%   daily - indicator daily/monthly data
%   drop_perc - drop characteristics with > drop_perc obs. missing
%              (default: 100 -- do not drop)
%   omit_prefixes - cell array of characteristics prefixes to be dropped
%   keeponly - keep only these characteristics
%       
% Output parameters:
%   dates - Dates
%   re - excess returns
  

function [dates, re, mkt, names, DATA] = load_managed_portfolios(filename, daily, drop_perc, omit_prefixes, keeponly)

    if nargin < 3
        drop_perc = 1;
    end
    if nargin < 4
        omit_prefixes = {};
    end
    if nargin < 5
        keeponly = {};
    end

    if daily 
        suffix = '_d';
        date_format = 'mm/dd/yyyy';
    else
        suffix = '';
        date_format = 'mm/yyyy';
    end
    
    DATA = readtable(filename); DATA.date = datenum2(DATA.date, date_format); 


    if ~isempty(keeponly)
        DATA = DATA(:, ['date' 'rme' keeponly]);
    else
        % drop characteristics with certain prefixes
        for i = 1:length(omit_prefixes)
            prefix = omit_prefixes{i};
            idx = strncmp(DATA.Properties.VariableNames, prefix, length(prefix))';
            DATA(:, idx) = [];
        end
    end
    
    % drop characteristics with > drop_perc percentage of obs missing
    dropidx = find(sum(isnan(DATA{:, 1:end}))'>size(DATA,1)*drop_perc);
    DATA(:, dropidx) = [];
    
    % drop dates with missing obs
    idx2keep = find(sum(isnan(DATA{:, 3:end}),2)==0);
    assert(length(idx2keep)>0.75*size(DATA,1), 'More than 25% of obs. need to be dropped!')
    DATA = DATA(idx2keep, :); 
    
    
    dates = DATA.date;
    mkt = DATA.rme;
    re = DATA{:, 4:end}; 
    names = DATA.Properties.VariableNames(4:end);

    % de-market
%     b = olsgmm(re, [ones(size(mkt,1),1) mkt], 0, -1);
%     re = re - mkt*b(2,:);    
end
