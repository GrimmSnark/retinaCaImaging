function [d,lwrCross,uprCross,lwrRef,uprRef, stateLevs] ...
    = transdursState(sig, polarityFlag, plotFlag, varargin)
%TRANSDURS Extract durations of bilevel waveform transitions
%   If polarityFlag is +1 return duration metrics for positive transitions
%   If polarityFlag is -1 return duration metrics for negative transitions
%   If polarityFlag is  0 return duration metrics for all transitions
%
%   If plotFlag is true, plot the locations of the upper and lower
%   crossings.
%
%   This function is for internal use only. It may be removed in the
%   future.
%
%   20230131 - Modified by M Savage to include stateLevel output

%   Copyright 2011-2021 The MathWorks, Inc.
%#codegen

% check to see if input is a row vector
needsTranspose = coder.internal.isConstTrue(isrow(sig));
if ~coder.internal.isConst(size(sig))
    coder.internal.assert(coder.internal.isConstTrue(isvector(sig)),...
        'signal:codegeneration:InputNotVector','x');
end
% extract leading numeric arguments
[x, t, n] = chktransargs(0, sig, varargin{:});

% extract trailing optional arguments
[tol, prl, stateLevs] = chktransopts(x, {'PctRefLevels'}, varargin{n:end});

% compute state boundaries and reference levels
% (set mid reference level to mean of upper and lower reference level)
meanPrl = (prl(1) + prl(2))/2;
[lwrBnd, uprBnd, lwrRef, midRef, uprRef] = ...
    pctreflev(stateLevs, tol, 100-tol, prl(1), meanPrl, prl(2));

% extract the transitions
tm = signal.internal.getTransitions(x, t, uprBnd, lwrBnd, uprRef, midRef, lwrRef);

% find the desired polarities
if polarityFlag<0
    idx = find(tm.Polarity<0);
    msgIdTag = 'FallTime';
elseif polarityFlag>0
    idx = find(tm.Polarity>0);
    msgIdTag = 'RiseTime';
else
    idx = 1:numel(tm.Polarity);
    msgIdTag = 'SlewRate';
end

dCol = tm.Duration(idx);
lwrCrossCol = tm.LowerCross(idx);
uprCrossCol = tm.UpperCross(idx);

if plotFlag
    coder.internal.assert(coder.target('MATLAB'),'signal:codegeneration:PlottingIsOnlyInMatlab')
    xdata = [lwrCrossCol lwrCrossCol uprCrossCol uprCrossCol]';
    ydata = repmat([lwrRef uprRef uprRef lwrRef]',1,numel(lwrCrossCol));
    [~, l] = signal.internal.plotInitialize( ...
        ['signal:internal:plotInitialize:PlotName' msgIdTag], x, t, ...
        xdata, ydata, ['signal:internal:plotInitialize:PlotPatch' msgIdTag]);
    l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendUpperCross', ...
        uprCrossCol, uprRef * ones(size(idx)),'rx ', ...
        'MarkerSize',10,'LineWidth',0.8);
    l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendLowerCross', ...
        lwrCrossCol, lwrRef * ones(size(idx)),'gx ', ...
        'MarkerSize',10,'LineWidth',0.8);
    signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [prl(1) NaN prl(2)]);
end

if needsTranspose
    d = dCol.';
    lwrCross = lwrCrossCol.';
    uprCross = uprCrossCol.';
else
    d = dCol;
    lwrCross = lwrCrossCol;
    uprCross = uprCrossCol;
end
