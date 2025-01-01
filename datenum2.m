function n = datenum2(arg1,arg2,arg3,h,min,s)

if isa(arg1, 'datetime')
    n = datenum(arg1.Year, arg1.Month, arg1.Day);
else
    switch nargin
        case 1
            n = datenum(arg1);
        case 2
            n = datenum(arg1,arg2);
        case 3
            n = datenum(arg1,arg2,arg3);
        otherwise
            n = datenum(arg1,arg2,arg3,h,min,s);
    end          
end