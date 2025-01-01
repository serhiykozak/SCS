function [names] = anomnames(anom)

    for i = 1:length(anom)
        a = anom{i};
        prefix = a(1:3);
        s = a(4:end);
        
        switch a(1:3)
            case 'rme'
                n = 'market';
            case 're_'
                n = s;
            case 'r2_'
                n = [s '$^2$'];
            case 'r3_'
                n = [s '$^3$'];
            case 'rX_'
                xsep = '__';
                idx = strfind(s, xsep);
                if isempty(idx)
                    xsep = '_'; % OLD convention
                    idx = strfind(s, xsep); 
                end
                n = [s(1:idx-1) '$\times$' s(idx+length(xsep):end)];
            otherwise
                switch prefix(1:2)
                    case 'r_'
                        n = a(3:end);
                    otherwise
                    n = a;
                end
        end
        
        n = strrep(n, 'I_DUM', 'ind');
%         idx = strfind(n, 'I_DUM');
%         if ~isempty(idx)
%             n = [n(1:idx-1) 'ind' n(idx+5:end)];
%         end
        
        n = strrep(n, '_', '\_');
        
        names{i} = lower(n);
    end

end
