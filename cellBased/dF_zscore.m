function [zScore, sIQR] = dF_zscore(dF)

% remove negative values
% dF_rezeroed = dF;
% dF_rezeroed( dF_rezeroed <0) = 0;
% percentile10 = prctile(dF_rezeroed, 10);
% percentile90 = prctile(dF_rezeroed, 90);
% 
% bottom10 = dF_rezeroed(dF_rezeroed<= percentile10 );
% top10 = dF_rezeroed(dF_rezeroed>= percentile90 );

baselinePercentile= prctile(dF, 20);
signalPercentile = prctile(dF, 99);

baseline = dF(dF<= baselinePercentile );
signal = dF(dF>= signalPercentile );

zScore = (mean(signal)- mean(baseline) )/ std(baseline);
sIQR = iqr(dF)/2;
end