function p = genpath_safe(d)
%genpath_safe genpath without slprj and SVN folders.
%   p = genpath_safe(d) is like the built-in genpath, except
%   that it attempts to remove slprj and SVN folders.

% Check input arguments
if isequal(nargin, 0)
    d = fullfile(matlabroot, 'toolbox'); % default
else
    d_exists = exist(d, 'dir') > 0;
    if ~d_exists
        id = 'genpath_safe:dirNotFound';
        msg = 'Directory not found.';
        ME = MException(id, msg);
        throwAsCaller(ME);
    end
end

% Use the built-in genpath to get a full subdirectory list.
full_path = genpath(d);
path_array = splitPathApart(full_path);

% Define the "illegal" patterns to look for at the start of a directory or
% file name. (Add to these as needed.)
illegal_patterns = {'.svn', '_svn', 'slprj', '.git'};
illegal_patterns = prepadWithFileSep(illegal_patterns);

% Remove the illegal patterns
fcn_handle = @(p)containsNoIllegalPatterns(p, illegal_patterns);
legal_idx = cellfun(fcn_handle, path_array);
path_array = path_array(legal_idx);

% Put everything back together again.
p = reassemblePath(path_array);

end


function p_array = splitPathApart(p)
pat = pathsep();
p_array = regexp(p, pat, 'split');
p_array = removeEmptyCells(p_array); % removes unwanted trailing cell
end


function s_out = removeEmptyCells(s_in)
anon_fcn = @(x)~isempty(x);
idx = cellfun(anon_fcn, s_in);
s_out = s_in(idx);
end


function patterns_out = prepadWithFileSep(illegal_patterns)
pat = filesep();
anon_fcn = @(x)horzcat(pat, x);
patterns_out = cellfun(anon_fcn, illegal_patterns, 'UniformOutput', false);
end


function success = containsNoIllegalPatterns(str, illegal_patterns)
all_success = cellfun(@(pat)doesNotContain(str, pat), illegal_patterns);
success = all(all_success);
end


function success = doesNotContain(str, illegal_pattern)
str_lower = lower(str);
pat_lower = lower(illegal_pattern);
locs = strfind(str_lower, pat_lower);
success = isempty(locs);
end


function p = reassemblePath(path_array)
pat = pathsep();
p = sprintf(['%s',pat],path_array{:});
p = p(1:end-numel(pat)); % removes unwanted trailing pathsep
end
