function [cfg, overridden, defaulted] = parse_config(cfg, CFG)
% Parse and process typical options in configuration string or structure
% Parameters:
%   cfg     - configuration to process
%   CFG     - defaults
% Returns:
%   cfg         - the parsed configuration
%   overridden  - fields that were explicitly set
%   defaulted   - fields that defaulted
%
%
% Copyright 2017, Serhiy Kozak. 
%
% First Version: 2017-07-21. Serhiy Kozak.
% This  Version: 2017-07-21. Serhiy Kozak.
%

% make sure cfg is struct even if empty
if isempty(cfg), cfg = struct; end;
    
% if provided string, parse it into a structure
if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;

% install defaults for whatever hasn't been specified up to now
[cfg, ~, overridden, defaulted] = check_params(cfg, CFG);


