% Demonstrates the Rea et al. 2005 model for melatonin suppression
function demo_rea2005model_monochromaticSpectra()

    % Petteri Teikari, 2013, LRC, Troy, NY, USA
    % teikap@rpi.edu
    close all    
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    
    %% SETTINGS
    
        % with 1 skips some tedious calculations
        loadFromMat = 0; 
    
        % FOLDERS
        path = initFolders();  
            
         % Styling
        style.fontName = 'latin modern roman';
        style.fontBaseSize = 9;
        style.fontLabelWeight = 'bold';

        % line plot colors, RGB-colors (divided by 255)
        style.colorGray = [0.2 0.2 0.2];
        style.colorPlot(1,:) = [0 0 0];
        style.colorPlot(2,:) = [0 1 1];
        style.colorPlot(3,:) = [0.87 0.49 0];
        style.colorPlot(4) = style.colorPlot(2);

        % FIGURE autosave options
        style.imgOutRes       = '-r300'; % dpi
        style.imgOutAntiAlias = '-a1';   % '-a0' least, '-a4' maximum
        style.imgOutautoSavePlot = 0;
        
        % Plot switches
        plotON.spectra = 0;
        plotON.components = 1; 
        plotON.criterionComparison = 0;
        
    %% MODEL PARAMETERS
    
        % REA et al. PARAMETERS
        simRes = 1; % nm resolution for spectra
        criteriaVector = [150 225 300 375 450]; % 300 is the default value
        FWHM = 10; % of the monochromatic light spectra
        spectralCrossover = 'sharp'; 
        % spectralCrossover = 'smooth'; 
        modelType = 'original';  
        % modelType = 'origWithCones';
        a_cone = 0.0;
        
        % test the model with different intensities  (W/m^2)
        % the vector min, max, and the spacing (resolution)
        irradMin = -4;
        irradMax = 3;
        irradRes = 0.25;       
        photDensityCnst = 0;
        
        % Ocular Media parameters
        age = [25 65];
        offset = 0.111; % 0.111 default for van de Kraats and van Norren 2007
        lambda = (380:1:780)'; % have to be the same as for light sources
        
        
        
    %% IMPORT THE Spectral Sensitivities of Photoreceptors        
    
        % Corneal Sensitivities
        cd(path.photoreception)
        S_cornea = import_CornealSensitivities(path.templates);
        cd(path.mainCode)
        
        % Retinal Sensitivities
        cd(path.photoreception)
        S_retina = import_RetinalSensitivities(path.nomogram);
        cd(path.mainCode)
        
        % Melanopsin parameters for bistable computation
        Rpeak   = 480; % use the one from Enezi et al., 2011
        Mpeak   = 587; % Mure et al., 2009         
        [alphaR_rel, alphaM_rel] = defineBistableMelanopsin(Rpeak,Mpeak,path); % small wrapper           
        
    %% DEFINE the MONOCHROMATIC SPECTRA    
    
        lambdaSimSPD = (min(lambda) : simRes : max(lambda))';

        % create monochromatic light, irradiance weighed
        [SPD.raw{1}, SPD.orig{1}, irrad] = createMonochromaticLightSPDs(lambda, lambdaSimSPD, FWHM, photDensityCnst, irradMin, irradMax, irradRes, path);

        % correct for ocular media transmission
        [SPD.filt, lensFilter] = correctForOcularMedia(lambdaSimSPD, SPD, irrad, age, lambda, offset, path);             
                    
    %% Compute responses and parameters
    
        tic
        if loadFromMat == 1
            load('MatSpectralSimul.mat')              
        else
            rea2005 = computeReaParameters(lambdaSimSPD, SPD.filt, irrad, criteriaVector, S_cornea, S_retina, spectralCrossover, modelType, a_cone, path);
            % save('MatSpectralSimul.mat', 'lambdaSimSPD', 'SPD', 'lensFilter', 'rea2005') % SAVE THE RESULTS to MAT
    
        end
        reaTime = toc
            
    %% PLOT The Results
    
    
        % for plots with only one criterion response, we will pick the
        % default one, quick'n'dirty check
        criterionIndex = find(criteriaVector == 300);
        if isempty(criterionIndex)
            error('You never used the default criterion value for calculating the responses')
        end
        
    
        %% For the spectral sensitivity
        if plotON.spectra == 1        
        
            % Linear
            linLog = 'lin';
            fig1 = plotJustTheSpectra(scrsz, style, lambdaSimSPD, SPD.filt, rea2005, criterionIndex, linLog);

                % autosave the figure      
                if style.imgOutautoSavePlot == 1
                    fileNameOut = ['reaModel2005_monochromaticSpectra_LIN_FWHM', num2str(FWHM), 'nm_res', num2str(simRes), 'nm_', spectralCrossover, '_', modelType,  '.png'];
                    try
                        export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                    catch err
                        warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                    end
                end

            % LOG
            linLog = 'log';
            fig2 = plotJustTheSpectra(scrsz, style, lambdaSimSPD, SPD.filt, rea2005, criterionIndex, linLog);

                % autosave the figure      
                if style.imgOutautoSavePlot == 1
                    fileNameOut = ['reaModel2005_monochromaticSpectra_LOG_FWHM', num2str(FWHM), 'nm_res', num2str(simRes), 'nm_', spectralCrossover, '_', modelType,  '.png'];
                    try
                        export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                    catch err
                        warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                    end
                end
                
        end
            
        %% For the components
        if plotON.components == 1
        
            % Linear
            linLog = 'lin';
            fig3 = plotWithComponentsMain(scrsz, style, lambdaSimSPD, SPD.filt, rea2005, criterionIndex, linLog);

                % autosave the figure      
                if style.imgOutautoSavePlot == 1
                    fileNameOut = ['reaModel2005_monochromaticComponents_', num2str(simRes), 'nm_irradRes', num2str(irradRes), '.png'];
                    try
                        export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                    catch err
                        warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                    end
                end

            % Logaritmic
            linLog = 'log';
            fig4 = plotWithComponentsMain(scrsz, style, lambdaSimSPD, SPD.filt, rea2005, criterionIndex, linLog);

                % autosave the figure      
                if style.imgOutautoSavePlot == 1
                    fileNameOut = ['reaModel2005_monochromaticComponents_LOG', num2str(simRes), 'nm_irradRes', num2str(irradRes), '.png'];
                    try
                        export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                    catch err
                        warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                    end                        
                end          
                
        end
        
        %% For the different criterion response values
        if plotON.criterionComparison == 1
            
            % Linear
            linLog = 'lin';
            fig5 = plotWithDifferentCriterionVectors(scrsz, style, lambdaSimSPD, SPD.filt, rea2005, criteriaVector, linLog);

                % autosave the figure      
                if style.imgOutautoSavePlot == 1
                    fileNameOut = ['reaModel2005_criterionResponses_', num2str(simRes), 'nm_irradRes', num2str(irradRes), '.png'];
                    try
                        export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                    catch err
                        warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                    end
                end
            
            
        end