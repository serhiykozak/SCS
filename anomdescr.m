function [descriptions] = anomdescr(anom)
    assert(iscell(anom) && ~isempty(anom), 'The first parameter must be a non-empty cell array!')

    % params
    n = length(anom);
    descriptions = cell(n, 1);
    
    % map from names to descriptions
    desc = characteristics_names_map();
    
    % translate all names
    for i = 1:length(anom)
        a = anom{i};
        prefix = a(1:3);
        s = a(4:end);
        
        switch a(1:3)
            case 'rme'
                n = 'Market';
            case 're_'
                n = read_desc(desc, s);
            case 'r2_'
                n = [read_desc(desc, s) '$^2$'];
            case 'r3_'
                n = [read_desc(desc, s) '$^3$'];
            case 'rX_'
                xsep = '__';
                idx = strfind(s, xsep);
                if isempty(idx)
                    xsep = '_'; % OLD convention
                    idx = strfind(s, xsep); 
                end
                n = [read_desc(desc, s(1:idx-1)) '$\times$' read_desc(desc, s(idx+length(xsep):end))];
            otherwise
                switch prefix(1:2)
                    case 'r_'
                        s = a(3:end);
                    otherwise
                        s = a;
                end
                
                n = read_desc(desc, s);
        end
        
        n = strrep(n, '_', '\_');
        
        descriptions{i} = n;
    end

end

% read description from the map, but do not throw an error if no desc is
% available
function [desc] = read_desc(map, s)
    if map.isKey(s)
        desc = map(s);
    else
        warning(['No description available for [' s ']'])
        desc = s;
    end
end

