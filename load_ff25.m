function [dates, ret, mkt, DATA, labels] = load_ff25(datapath, daily, t0, tN)

    if nargin < 3
        t0 = 0;
    end
    if nargin < 4
        tN = inf;
    end

    if daily 
        ffact5 = 'F-F_Research_Data_Factors_daily.csv';
        ff25 = '25_Portfolios_5x5_Daily_average_value_weighted_returns_daily.csv';
        %fmom = 'F-F_Momentum_Factor_daily.csv';
    else
        ffact5 = 'F-F_Research_Data_Factors.csv';
        ff25 = '25_Portfolios_5x5_average_value_weighted_returns_monthly.csv';
        %fmom = 'F-F_Momentum_Factor.csv';
    end
    
    def_date_format = 'yyyy/mm/dd';

    DATA = readtable([datapath ffact5]); DATA.Date = datenum(DATA.Date, def_date_format); 
    DATA(DATA.Date<t0 | DATA.Date>tN, :) = [];
    
    RET = readtable([datapath ff25]); RET.Date = datenum(RET.Date, def_date_format); 
    
    DATA = innerjoin(DATA, RET, 'keys', 1); 

    dates = DATA.Date;
    mkt = DATA.Mkt_RF/100;
    ret = DATA{:, 6:30}/100 - DATA.RF/100;
    labels = RET.Properties.VariableNames(2:end);
end


%     'date'
%     'Mkt_RF'
%     'SMB'
%     'HML'
%     'RF'
%     'SMALLLoBM_DATA'
%     'ME1BM2_DATA'
%     'ME1BM3_DATA'
%     'ME1BM4_DATA'
%     'SMALLHiBM_DATA'
%     'ME2BM1_DATA'
%     'ME2BM2_DATA'
%     'ME2BM3_DATA'
%     'ME2BM4_DATA'
%     'ME2BM5_DATA'
%     'ME3BM1_DATA'
%     'ME3BM2_DATA'
%     'ME3BM3_DATA'
%     'ME3BM4_DATA'
%     'ME3BM5_DATA'
%     'ME4BM1_DATA'
%     'ME4BM2_DATA'
%     'ME4BM3_DATA'
%     'ME4BM4_DATA'
%     'ME4BM5_DATA'
%     'BIGLoBM_DATA'
%     'ME5BM2_DATA'
%     'ME5BM3_DATA'
%     'ME5BM4_DATA'
%     'BIGHiBM_DATA'
%     'SMALLLoBM_AFSIZE'
%     'ME1BM2_AFSIZE'
%     'ME1BM3_AFSIZE'
%     'ME1BM4_AFSIZE'
%     'SMALLHiBM_AFSIZE'
%     'ME2BM1_AFSIZE'
%     'ME2BM2_AFSIZE'
%     'ME2BM3_AFSIZE'
%     'ME2BM4_AFSIZE'
%     'ME2BM5_AFSIZE'
%     'ME3BM1_AFSIZE'
%     'ME3BM2_AFSIZE'
%     'ME3BM3_AFSIZE'
%     'ME3BM4_AFSIZE'
%     'ME3BM5_AFSIZE'
%     'ME4BM1_AFSIZE'
%     'ME4BM2_AFSIZE'
%     'ME4BM3_AFSIZE'
%     'ME4BM4_AFSIZE'
%     'ME4BM5_AFSIZE'
%     'BIGLoBM_AFSIZE'
%     'ME5BM2_AFSIZE'
%     'ME5BM3_AFSIZE'
%     'ME5BM4_AFSIZE'
%     'BIGHiBM_AFSIZE'
%     'SMALLLoBM_FNUM'
%     'ME1BM2_FNUM'
%     'ME1BM3_FNUM'
%     'ME1BM4_FNUM'
%     'SMALLHiBM_FNUM'
%     'ME2BM1_FNUM'
%     'ME2BM2_FNUM'
%     'ME2BM3_FNUM'
%     'ME2BM4_FNUM'
%     'ME2BM5_FNUM'
%     'ME3BM1_FNUM'
%     'ME3BM2_FNUM'
%     'ME3BM3_FNUM'
%     'ME3BM4_FNUM'
%     'ME3BM5_FNUM'
%     'ME4BM1_FNUM'
%     'ME4BM2_FNUM'
%     'ME4BM3_FNUM'
%     'ME4BM4_FNUM'
%     'ME4BM5_FNUM'
%     'BIGLoBM_FNUM'
%     'ME5BM2_FNUM'
%     'ME5BM3_FNUM'
%     'ME5BM4_FNUM'
%     'BIGHiBM_FNUM'
%     'year'
%     'Date'
%     'SMALLLoBM_SUMBE'
%     'ME1BM2_SUMBE'
%     'ME1BM3_SUMBE'
%     'ME1BM4_SUMBE'
%     'SMALLHiBM_SUMBE'
%     'ME2BM1_SUMBE'
%     'ME2BM2_SUMBE'
%     'ME2BM3_SUMBE'
%     'ME2BM4_SUMBE'
%     'ME2BM5_SUMBE'
%     'ME3BM1_SUMBE'
%     'ME3BM2_SUMBE'
%     'ME3BM3_SUMBE'
%     'ME3BM4_SUMBE'
%     'ME3BM5_SUMBE'
%     'ME4BM1_SUMBE'
%     'ME4BM2_SUMBE'
%     'ME4BM3_SUMBE'
%     'ME4BM4_SUMBE'
%     'ME4BM5_SUMBE'
%     'BIGLoBM_SUMBE'
%     'ME5BM2_SUMBE'
%     'ME5BM3_SUMBE'
%     'ME5BM4_SUMBE'
%     'BIGHiBM_SUMBE'




