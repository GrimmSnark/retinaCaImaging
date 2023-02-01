function varargout = chktransopts(x, opts, varargin)
%CHKTRANSOPTS check arguments for bilevel waveform transitions
%
%   out position   option                     default
%       1          'Tolerance'                2
%       2          'MidPct' or 'PctRefLevels' 50 or [10 90]
%       3          'StateLevels'              statelevels(x)
%       4          'Polarity'                 'positive'
%    4 and 5       'Aberration'               'postshoot' and 3
%
%   Duplicated/unrecognized options generate an error
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2021 The MathWorks, Inc.
%#codegen

if ~coder.target('MATLAB')
    for i = 1:length(varargin)
        if ischar(varargin{i}) || isStringScalar(varargin{i})
            coder.internal.assert(coder.internal.isConst(varargin{i}),...
                'signal:codegeneration:OptionInputsMustBeConstant');
        end
    end
end

% idxArray : boolean vector having length equal to the number of elements in
% varargin. All elements of idxArray are initially true. After the parsing
% is finished, if idxArray(i) remains true, then varargin{i} was an illegal
% input.

idxArray = true(1,length(varargin));
[tol,idxArray] = extractTol(varargin,idxArray);
varargout{1} = tol;

if any(strcmpi(opts,'MidPct'))
    [varargout{2}, idxArray] = extractMprl(tol,varargin,idxArray);
else
    [varargout{2}, idxArray] = extractPrl(tol,varargin,idxArray);
end

[varargout{3}, idxArray] = extractStateLevs(x,varargin,idxArray);

if any(strcmpi(opts,'Polarity'))
    [varargout{4}, idxArray] = extractPol(varargin,idxArray);
end

if any(strcmpi(opts,'Region'))
    [varargout{4}, idxArray] = extractRgn(varargin,idxArray);
    [varargout{5}, idxArray] = extractFactor(varargin,idxArray);
end

for i = 1:length(idxArray)
    coder.internal.errorIf(idxArray(i) && (ischar(varargin{i}) || isStringScalar(varargin{i})),....
        'signal:chktransopts:UnexpectedOption',varargin{i});
    coder.internal.errorIf(idxArray(i),...
        'signal:chktransopts:UnexpectedInput')
end
end

function [val,idxArray,found] = findprop(prop,default,idxArray,varargin)

found = false;
coder.unroll();
for i = 1:length(varargin)
    if strcmpi(varargin{i},prop)
        coder.internal.errorIf(found,'signal:chktransopts:DuplicateProperty',prop);
        coder.internal.errorIf(isempty({varargin{i+1:end}}),'signal:chktransopts:MissingValue',prop);
        found = true;
        val = varargin{i+1};
        idxArray(i) = false;
        idxArray(i+1) = false;
    end
end

if ~found
    val = default;
end

end

function [tol,idxArray] = extractTol(p,idxArray)
[tol, idxArray] = findprop('tolerance',2,idxArray,p{:});
validateattributes(tol, ...
    {'double'},{'real','finite','positive','scalar','<',50},'','TOL');
tol = tol(1);
end

function [mprl,idxArray] = extractMprl(tol, p,idxArray)
[mprl, idxArray,foundmprl] = findprop('midpercentreferencelevel',50,idxArray,p{:});
if ~foundmprl
    [mprl, idxArray] = findprop('midpct',50,idxArray,p{:});
end
validateattributes(mprl, ...
    {'double'},{'real','finite','scalar'},'','MIDREFERENCELEVEL');

mprl = mprl(1);
coder.internal.assert(tol < mprl,'signal:chktransopts:TolMustBeLessThanMidPctRefLev');
coder.internal.assert(mprl+tol < 100,'signal:chktransopts:TolAndMidPctRefLevTooBig');
end

function [prl,idxArray] = extractPrl(tol, p,idxArray)
[prl, idxArray,foundprl] = findprop('percentreferencelevels',[10 90],idxArray,p{:});
if ~foundprl
    [prl, idxArray] = findprop('pctreflevels',[10 90],idxArray,p{:});
end
validateattributes(prl,{'double'},{'real','finite','size',[1 2]}...
    ,'','PCTREFLEVELS');

coder.internal.assert(prl(1) < prl(2),'signal:chktransopts:LowerMustBeLessThanUpperPctRefLev')
coder.internal.assert(tol < prl(1),'signal:chktransopts:TolMustBeLessThanLwrPctRefLev');
coder.internal.assert(prl(2)+tol <100,'signal:chktransopts:TolAndUprPctRefLevTooBig');
end

function [stateLevs, idxArray] = extractStateLevs(x,p,idxArray)
[stateLevs, idxArray,foundLevel] = findprop('statelevels',[0 1],idxArray,p{:});
if foundLevel
    validateattributes(stateLevs,{'double'},{'real','finite','size',[1 2]},...
        '','STATELEVELS');
    coder.internal.assert(stateLevs(1) < stateLevs(2),...
        'signal:chktransopts:LowerMustBeLessThanUpperStateLev')
else
    stateLevs = statelevels(x);
end
end

function [f, idxArray] = extractFactor(p,idxArray)
[f, idxArray] = findprop('seekfactor',3,idxArray,p{:});
validateattributes(f, ...
    {'double'},{'real','finite','scalar','positive'},'','FACTOR');
f = f(1);
end

function [polFlag, idxArray] = extractPol(p,idxArray)
[polStr, idxArray] = findprop('polarity','positive',idxArray,p{:});

coder.internal.assert(strcmpi(polStr,'positive') || strcmpi(polStr,'negative'),...
    'signal:chktransopts:IllegalBooleanFlagValue','POLARITY','Positive','Negative')

if strcmpi(polStr,'positive')
    polFlag = 1;
else
    polFlag = -1;
end
end

function [rgn, idxArray] = extractRgn(p,idxArray)
[rgn, idxArray] = findprop('region','postshoot',idxArray,p{:});
coder.internal.assert(strcmpi(rgn,'preshoot') || strcmpi(rgn,'postshoot'),...
    'signal:chktransopts:IllegalBooleanFlagValue','REGION','Preshoot','Postshoot')
end
