function fayeCode(exStructPath)

%% open FIJI
% initalize MIJI and get ROI manager open
intializeMIJ;
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% load metaStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;


%% analyse each cell

for c= 1:exStruct.cellCount
    dFTrace = exStruct.cells.dF(c,:);

    % min peak height
    dF_SD = std(dFTrace);
    dFMean = mean(dFTrace);
    minHeight = dFMean + (dF_SD *2);
    minProm = 2.5 * dF_SD;

    % spike detection
       findpeaks(dFTrace, "MinPeakProminence",minProm);
%     [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace,"MinPeakHeight",minHeight, "MinPeakProminence",minProm);
        % get half width with peaks

        % get tau

 end

end