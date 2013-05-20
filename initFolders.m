function path = initFolders()

    mFilename = 'initFolders'; 
    pathF = mfilename('fullpath');
    path.mainCode = strrep(pathF, mFilename, '');
    cd(path.mainCode); cd ..; cd ..;
    path.lightLab = pwd;
    
    path.library        = path.mainCode;
    
%     path.nomogram       = fullfile(path.lightLab, 'lib', 'nomogram');
%     path.lightSources   = fullfile(path.lightLab, 'database', 'LightSources', 'naturalSources');
%     path.bistability    = fullfile(path.lightLab, 'lib', 'bistability');            
%     path.nomogram       = fullfile(path.lightLab, 'lib', 'nomogram');   
%     path.colorimetry    = fullfile(path.lightLab, 'lib', 'colorimetry');   
%     path.photoreception = fullfile(path.lightLab, 'lib', 'photoreception');
%     path.common         = fullfile(path.lightLab, 'lib', 'common');   
%     path.templates      = fullfile(path.lightLab, 'database', 'Templates');
%     path.ocularMedia    = fullfile(path.lightLab, 'lib', 'ocularmedia');
      
    
    % "Original CS files" from Andy
    path.rea2005orig    = fullfile(path.mainCode, 'Rea2005');
    
    % Data folders
    path.dataOut        = fullfile(path.mainCode, 'dataOut');
    path.dataIn         = fullfile(path.mainCode, 'dataIn'); 
    path.figuresOut     = fullfile(path.mainCode, 'figuresOut');
    
    % Subfunctions needed for some of the calculations  
    path.photoreception = fullfile(path.mainCode, 'photoreception');
    path.ocularMedia    = fullfile(path.mainCode, 'ocularMedia');
    path.nomogram       = fullfile(path.mainCode, 'nomogram');
    path.templates      = fullfile(path.mainCode, 'templates');

    
    
    cd(path.mainCode)