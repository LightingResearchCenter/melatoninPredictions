function plot_sensitivityPlanes()

    %% Parameters
    
        % FOLDERS
            mFilename = 'plot_sensitivityPlanes'; 
            pathF = mfilename('fullpath');
            path.mainCode = strrep(pathF, mFilename, '');
            cd(path.mainCode); cd ..; cd ..; cd ..;
            path.lightLab = pwd;
            path.library        = fullfile(path.lightLab, 'lib');
            path.nomogram       = fullfile(path.lightLab, 'lib', 'nomogram');
            path.lightSources   = fullfile(path.lightLab, 'database', 'LightSources', 'naturalSources');
            path.bistability    = fullfile(path.lightLab, 'lib', 'bistability');            
            path.nomogram       = fullfile(path.lightLab, 'lib', 'nomogram');   
            path.colorimetry    = fullfile(path.lightLab, 'lib', 'colorimetry');   
            path.photoreception = fullfile(path.lightLab, 'lib', 'photoreception');
            path.common         = fullfile(path.lightLab, 'lib', 'common');   
            path.templates      = fullfile(path.lightLab, 'database', 'Templates');
            path.ocularMedia    = fullfile(path.lightLab, 'lib', 'ocularmedia');
            path.figuresOut     = fullfile(path.lightLab, 'figures');
            
        % wavelength vector
        xRes = 1;
        xLimits = [380 780];
        lambda = (xLimits(1):xRes:xLimits(2))';
        FWHM = 1;        
            
        % Go through the possible monochromatic light stimuli
        for i = 1 : length(lambda)
            
            SPD{i} = monochromaticLightAsGaussian(lambda(i), FWHM, xRes, xLimits)
        
            for j = 1 : length(irradiances)
                
                % CALCULATE THE CLA
                [rea2005.CLA{i}{j}, rea2005.CLAComp{i}{j}, rea2005.CLASpec{i}{j}, oppFlag(i,j)] = ...
                                        SPD{i}, S_cornea, S_retina, lambMin, lambMax, lambRes);
                                    
                % CALCULATE THE CS
                rea2005.CS{i}{j} = CSCalc(rea2005.CLA{i}{j}); 
                
            end
            
        end