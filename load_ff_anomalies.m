function [dates ret mkt DATA] = load_ff_anomalies(datapath, daily, t0, tN)

    if nargin < 3
        t0 = 0;
    end
    if nargin < 4
        tN = inf;
    end

    if daily 
%         suffix = 'd';
        ffact5 = 'F-F_Research_Data_5_Factors_2x3_daily.csv';
        fmom = 'F-F_Momentum_Factor_daily.csv';
    else
%         suffix = 'm';
        ffact5 = 'F-F_Research_Data_5_Factors_2x3.csv';
        fmom = 'F-F_Momentum_Factor.csv';
    end
    
    def_date_format = 'yyyy/mm/dd';

%     DATA = readtable([datapath 'fact5_' suffix '.csv']); DATA.date = datenum(DATA.date, def_date_format); 
    DATA = readtable([datapath ffact5]); DATA.date = datenum(DATA.Date, def_date_format); 
    DATA(find(DATA.date<t0 | DATA.date>tN), :) = [];

%     MOM = readtable([datapath 'mom_' suffix '.csv'], 'Delimiter', {','}); MOM.date = datenum(MOM.date, def_date_format); 
    MOM = readtable([datapath fmom], 'Delimiter', {','}); MOM.date = datenum(MOM.Date, def_date_format); 
    
    DATA = innerjoin(DATA, MOM, 'keys', 'date'); 


    dates = DATA.date;
    ret = DATA{:, {'SMB', 'HML', 'Mom', 'RMW' ,'CMA'}}/100; 
    mkt = DATA.Mkt_RF/100;
    
    % de-market
%     b = olsgmm(ret, [ones(size(mkt,1),1) mkt], 0, -1);
%     ret = ret - mkt*b(2,:);    
end

