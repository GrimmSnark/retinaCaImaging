function [fittedY, yBaselined, coefficients] =fitExpCurveGPU(xData, yData)
% Fits an expoential decay curve on the data on matlab gpu
%
% Input- xData: x data in gpu Array
%        yData: y data in gpu array (pixel x timepoint array)
%
% Output- fittedY - the fitted exp curve y data (GPU array)
%         yBaselined - the baseline subtracted y data (GPU array)
%         coefficients - the estimated parameters for the fit

%%
% Take the logarithm of the y-values

logY = log(yData);
% logY = log(complex(yData));

% remove any inf values
logY(logY == -inf) = 0;

% Fit a linear polynomial to the transformed data
coefficients = polyfit(xData, logY, 1);

% Extract the slope and intercept from the coefficients
slope = coefficients(1);
intercept = coefficients(2);

% Compute the fitted exponential decay curve
fittedY = exp(intercept) * exp(slope * xData);
fittedY = fittedY* 0.95; % scale down to allow for postive trace


% fittedY = gather(fittedY);

yBaselined = yData - fittedY;

% add offset to avoid weird artefacts
if min(yBaselined) < 1000
    yBaselined = yBaselined+ abs(min(yBaselined)) + 1000;
end

% 
% % Plot the original data and the fitted curve
% figure;
% plot(x, y, 'ro', 'MarkerSize', 8); % Original data
% hold on;
% plot(x, fittedY, 'b-', 'LineWidth', 2); % Fitted curve
% plot(x,yBaselined2, 'LineWidth', 2);
% xlabel('x');
% ylabel('y');
% legend('Original Data', 'Fitted Curve');
% title('Exponential Decay Curve Fitting');

end