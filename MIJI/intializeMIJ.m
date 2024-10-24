% sets paths for miji information: will need to be personalised for your
% computer

% Set plot defaults
set(0,'defaultAxesTickDir','out')
set(0,'defaultAxesTickDirMode','manual')
set(0,'defaultAxesBox','off')
set(0,'DefaultFigureWindowStyle','normal') %docked or normal
set(0,'DefaultFigureColor','w')
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');

try  ij.IJ.getInstance().toFront();
    
    
catch ME
    
    % Add ImageJ to working directory
     javaaddpath([matlabroot '\java\jar\mij.jar']);
        javaaddpath([matlabroot '\java\jar\ij-1.54f.jar']);
%       javaaddpath([matlabroot '\java\jar\ij.jar']);
    
    % Add ImageJ plugins to the current path
    fijiPath = 'E:\PostDoc_Docs\Fiji.app\';
    javaaddpath([fijiPath '\plugins'])
    javaaddpath([fijiPath '\macros'])
%     javaaddpath([fijiPath 'plugins\BIJ_\bij.jar'])
    javaaddpath([fijiPath 'plugins\Cell_Magic_Wand_Tool.jar'])
    javaaddpath([fijiPath 'jars\microglia-morphometry-0.5.3.jar']);
    javaaddpath([fijiPath 'jars\commons-pool2-2.11.1.jar']);
    javaaddpath([fijiPath 'jars\labkit-ui-0.3.11.jar']);
    javaaddpath([fijiPath 'jars\labkit-pixel-classification-0.1.17.jar']);

%     javaaddpath([fijiPath 'plugins\Image_Stabilizer\'])
    % javaaddpath([fijiPath 'plugins\bUnwarpJ_.jar']);
    addpath([fijiPath 'scripts']);
    
    % Startup ImageJ
    ImageJ;
    
    MIJ.run('Install...', ['install=[' fijiPath '/macros/StartupMacros.fiji.ijm]']);
    
end

clear currentJavaPath
% IJ =ij.IJ();