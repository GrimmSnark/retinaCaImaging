function [fittedExpCurve, yBaselined, estimates, model] = fitexpCurve(xdata, ydata)
% Call fminsearch with a random starting point.
% start_point = rand(1, 3);
% guessed the following starting point;
start_point = [10 10 0 0.01];
model = @expfun;

estimates = fminsearch(model, start_point, optimset('MaxFunEvals',1000000, 'MaxIter', 1000000));


%% create exp baseline curve

A=estimates(1);
B=estimates(2);
lambda=estimates(3);
C=estimates(4);
fittedExpCurve = B + A .* exp(-lambda * xdata)+C*xdata;
fittedExpCurve = fittedExpCurve* 0.97; % scale down to allow for postive trace

%% create corrected baseline

yBaselined = ydata - fittedExpCurve;

% add offset to avoid weird artefacts
if min(yBaselined) < 1000
    yBaselined = yBaselined+ abs(min(yBaselined)) + 1000;
end

%% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurvePlusB] = expfun(params)
        A = params(1);
        B = params(2);
        lambda = params(3);
        C= params(4);
        FittedCurvePlusB = B + A .* exp(-lambda * xdata)+ C.*xdata;
        ErrorVector = FittedCurvePlusB - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end