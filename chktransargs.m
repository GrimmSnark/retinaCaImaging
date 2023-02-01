function [xOut, t, n, d] = chktransargs(needDelay, x, varargin)
%CHKTRANSARGS check arguments for bilevel waveform transitions
%
%   Validates and extracts (numeric) inputs of the form:
%       (X,     <optional string arguments>)
%       (X, FS, <optional string arguments>)
%       (X, T,  <optional string arguments>)
%
%   If the needDelay argument is true then extract inputs of the form:
%       (X, D,     <optional string arguments>)
%       (X, FS, D, <optional string arguments>)
%       (X, T, D,  <optional string arguments>)
%
%   It will return the sample values, x, and sample instants, t,
%   that correspond to the input vector as well as the index, n, into
%   varargin that specifies the first [unhandled] string argument.
%   It will also return the delay parameter, d, if it is requested.
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2021 The MathWorks, Inc.
%#codegen

    validateattributes(x,{'double'},{'real','finite','vector'}, '','X');
    coder.internal.assert(numel(x) >= 2,...
                          'signal:chktransargs:MustBeMultiElementVector','X');
    xOut = x(:);

    n = getNumExtraArgs(varargin{:});

    coder.internal.assert(n <= 1+needDelay,...
                          'signal:chktransargs:TooManyNumericArgs',2+needDelay)
    if needDelay
        coder.internal.assert(n~=0,'signal:chktransargs:MissingParameter','D');
        validateattributes(varargin{n},{'double'},...
                           {'real','finite','positive','scalar'},'','D');
        d = varargin{n}(1);
    end
    t = getTimeVector(x, varargin{1:n-needDelay});
    n = n+1;
end

function n = getNumExtraArgs(varargin)
    coder.inline('always');
    % get the number of extra numeric arguments for transition measurements.
    n = 0;
    coder.unroll();
    for i = 1:nargin
        if isnumeric(varargin{i})
            n = n + 1;
        else
            break;
        end
    end
end

function t = getTimeVector(x, varargin)
    coder.inline('always');

    if nargin==1
        t = (1:numel(x))';
    elseif isscalar(varargin{1})
        fs = varargin{1};
        validateattributes(fs,{'double'},{'real','finite','positive'},'','FS');
        t = (0:numel(x)-1)' / fs(1);
    else
        validateattributes(varargin{1},{'double'},...
                           {'real','finite','vector','increasing'},'','T');
        t = reshape(varargin{1},[],1);
        coder.internal.assert(numel(x) == numel(t),'signal:chktransargs:LengthMismatch')
    end
end
