function c = nanconvn(a, k, varargin)
%Full credit to Benjamin Kraus for making this code. This implementation in N dimensions simply
%extends the Matlab Exchange entry https://www.mathworks.com/matlabcentral/fileexchange/41961-nanconv


% NANCONVN Convolution in ND ignoring NaNs.
%   C = NANCONVN(A, K) convolves A and K, correcting for any NaN values
%   in the input vector A. The result is the same size as A (as though you
%   called 'conv' or 'conv2' with the 'same' shape).
%
%   C = NANCONVN(A, K, 'param1', 'param2', ...) specifies one or more of the following:
%     'edge'     - Apply edge correction to the output.
%     'noedge'   - Do not apply edge correction to the output (default).
%     'nanout'   - The result C should have NaNs in the same places as A.
%     'nonanout' - The result C should have ignored NaNs removed (default).
%                  Even with this option, C will have NaN values where the
%                  number of consecutive NaNs is too large to ignore.
%
%   NANCONV works by running 'convn' either two or three times. The first
%   time is run on the original input signals A and K, except all the
%   NaN values in A are replaced with zeros. The 'same' input argument is
%   used so the output is the same size as A. The second convolution is
%   done between a matrix the same size as A, except with zeros wherever
%   there is a NaN value in A, and ones everywhere else. The output from
%   the first convolution is normalized by the output from the second 
%   convolution. This corrects for missing (NaN) values in A, but it has
%   the side effect of correcting for edge effects due to the assumption of
%   zero padding during convolution. When the optional 'noedge' parameter
%   is included, the convolution is run a third time, this time on a matrix
%   of all ones the same size as A. The output from this third convolution
%   is used to restore the edge effects. The 'noedge' parameter is enabled
%   by default so that the output from 'nanconv' is identical to the output
%   from 'conv2' when the input argument A has no NaN values.
%
% See also conv, conv2, convn
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2013, Benjamin Kraus

% Added feature of ND convolution by Fernando Zigunov (zigunov.com).
% Copyright (c) 2020, Fernando Zigunov


% Process input arguments
for arg = 1:nargin-2
    switch lower(varargin{arg})
        case 'edge'; edge = true; % Apply edge correction
        case 'noedge'; edge = false; % Do not apply edge correction
        case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
        case 'nanout'; nanout = true; % Include original NaNs in the output.
        case 'nonanout'; nanout = false; % Do not include NaNs in the output.
    end
end
% Apply default options when necessary.
if(exist('edge','var')~=1); edge = false; end
if(exist('nanout','var')~=1); nanout = false; end
if(exist('shape','var')~=1); shape = 'same';
elseif(~strcmp(shape,'same'))
    error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
end
% Get the size of 'a' for use later.
sza = size(a);

% Flat function for comparison.
o = ones(size(a));
% Flat function with NaNs for comparison.
on = ones(size(a));
% Find all the NaNs in the input.
n = isnan(a);
% Replace NaNs with zero, both in 'a' and 'on'.
a(n) = 0;
on(n) = 0;
% Check that the filter does not have NaNs.
if(any(isnan(k)));
    error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
end
% Calculate what a 'flat' function looks like after convolution.
if(any(n(:)) || edge)
    flat = convn(on,k,shape);
else flat = o;
end
% The line above will automatically include a correction for edge effects,
% so remove that correction if the user does not want it.
% JCWS commented:
% if(any(n(:)) && ~edge); flat = flat./conv2(o,k,shape); end
% JCWS new
if(any(n(:)) && ~edge); flat = flat./convn(o,k,shape); end

% Do the actual convolution
c = convn(a,k,shape)./flat;
% If requested, replace output values with NaNs corresponding to input.
if(nanout); c(n) = NaN; end
end
