function [s, flag, overridden, defaulted] = check_params(s, defaults, required)
%Verifies parameter structure and sets defaults for optional parameters
%  [S] = CHECKPARAMS(S, DEFAULTS, REQUIRED)
%  Verifies a parameter structure for required fields and sets defaults for
%  optional fields. If any required field is missing, this will be signaled
%  with -1 on the FLAG. If no FLAG is required on the output, an error will
%  be fired.
%  
%  Parameters:
%   S           - the parameters structure to check
%   DEFAULTS    - the default values for optional parameters, a structure
%               with the same field names for parameters as those expected
%               in S. Can contain nested structures. Can be an empty matrix,
%               in which case no defaults are installed but the required fields
%               are verified.
%   REQUIRED    - the names of the required parameters, cell array of
%               strings
%
%  Note that type verifications are not performed; e.g. if a numeric is in
%  place of a structure, it will be handled as a structure and as a result
%  an (unaddressed) error will occur.
%
%  Returns:
%   S           - the parameters structure with the missing values replaced
%               by defaults
%   FLAG        - signals whether initialization finished correctly or not
%   OVERRIDDEN  - fields that were explicitly set
%   DEFAULTED   - fields that defaulted
%
%
% Copyright 2017, Serhiy Kozak. 
%
% First Version: 2017-07-21. Serhiy Kozak.
% This  Version: 2017-07-21. Serhiy Kozak.
%


REQARGIN = 2;
if nargin < REQARGIN + 1,
    required = {};
end

overridden = {}; defaulted = {};

for i = 1 : length(required)
    if ~isfield(s, required{i})
        if nargout > 1
            flag = -1;
            return;
        else 
            error(['Field ' required{i} ' is required']);
        end
    end
end

flag = 1;       % required OK

if isempty(defaults), return; end

fld = fieldnames(defaults);
for i = 1 : length(fld)
    if ~isfield(s, fld{i})            % assign defaults
        s.(fld{i}) = defaults.(fld{i});
        defaulted{end+1} = fld{i};
    else
        overridden{end+1} = fld{i};
        % stop at 1st level structures
%         % process recursively for defaults
%         if isstruct(s.(fld{i})),
%           s.(fld{i}) = checkparams(s.(fld{i}), defaults.(fld{i}));
%         end
    end
end

% END checkparams RETURNING s, flag ====================================