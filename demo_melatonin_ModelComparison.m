% Demonstrates the Rea et al. 2005 model for melatonin suppression
function demo_melatonin_ModelComparison()

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
        style.imgOutautoSavePlot = 1;
        
    %% MODEL PARAMETERS
    
        % REA et al. PARAMETERS
        simRes = 1; % nm resolution for spectra
        criteriaVector = [300 600];
        FWHM = 10; % of the monochromatic light spectra
        spectralCrossover = 'sharp'; 
        % spectralCrossover = 'smooth'; 
        modelType = 'original';  
        modelType = 'origWithCones';
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
    
        if loadFromMat == 1            
            load(fullfile(path.dataIn, 'MatModelCompar.mat'))            
            
        else
            
            % Rea et al. 2005
            %load(fullfile(path.dataIn, 'MatSpectralSimul.mat')) 
            rea2005 = computeReaParameters(lambdaSimSPD, SPD.filt, irrad, criteriaVector, S_cornea, S_retina, spectralCrossover, modelType, a_cone, path);          
            
             % Takahashi et al. 2011
            cd(path.photoreception)
            [pupil, melSuppr] = compute_TakahashiCurve(lambdaSimSPD, SPD.filt, irrad, path.templates);
            cd(path.mainCode)
            
            % save(fullfile(path.dataIn, 'MatModelCompar.mat'), 'lambdaSimSPD', 'SPD', 'lensFilter', 'rea2005', 'pupil', 'melSuppr') % SAVE THE RESULTS to MAT
            
        end
            
        for i = 1 : length(lensFilter)

            % Gall et al. 2004 / Lang et al. 2011
            cd(path.photoreception)
            gallDIN{i} = create_GallDIN(lambda, lensFilter{i}, path.templates);

            % Enezi et al. 2011
            cd(path.photoreception)
            melanopic{i} = create_MelanopicCurve(lambda, lensFilter{i}, path.templates);               

        end
        cd(path.mainCode)  
           
    
        
            
    %% PLOT The Results
    
            
        normStyle = {'toUnity'; 'to25years'; 'toLongWavelength'};
        linLog    = {'lin'}; %; 'log'};
        yStr = {'Rea/LRC'; 'Takahashi et al. 2011'; 'Gall/DIN, C(\lambda)'; 'Melanopic, M\phi(\lambda)'};
        titleStr = {'Rea/LRC'; 'Takahashi et al. 2011'; 'Gall/DIN, C(\lambda)'; 'Melanopic, M\phi(\lambda)'};
        handles = [];
        
        for i = 1 : length(normStyle)
        
            for k = 1 : length(linLog)
                
                fig = figure('Color','w');
                    rows = 2;
                    cols = 2;


                    %% REA MODEL
                    j = 1;
                    sp(j) = subplot(rows,cols,j);
                          [brainardStats{j}, thapanStats{j}, najjarStats{j}] = plotModelComparison(j, fig, scrsz, style, lambdaSimSPD, rea2005.ICLAeffNorm, normStyle{i}, linLog{k}, yStr{j}, titleStr{j}, path, handles);
                          
                    %% TAKAHASHI
                    j = 2;
                    sp(j) = subplot(rows,cols,j);
                        [brainardStats{j}, thapanStats{j}, najjarStats{j}] = plotModelComparison(j, fig, scrsz, style, lambda, melSuppr, normStyle{i}, linLog{k}, yStr{j}, titleStr{j}, path, handles);

                    %% GALL / DIN curve, C(lambda)
                    j = 3;
                    sp(j) = subplot(rows,cols,j);
                        [brainardStats{j}, thapanStats{j}, najjarStats{j}] = plotModelComparison(j, fig, scrsz, style, lambda, gallDIN, normStyle{i}, linLog{k}, yStr{j}, titleStr{j}, path, handles);

                    %% Melanopic, M_omega(lambda)
                    j = 4;
                    sp(j) = subplot(rows,cols,j);
                        [brainardStats{j}, thapanStats{j}, najjarStats{j}] = plotModelComparison(j, fig, scrsz, style, lambda, melanopic, normStyle{i}, linLog{k}, yStr{j}, titleStr{j}, path, handles);
            
                        % autosave the figure      
                        if style.imgOutautoSavePlot == 1
                            fileNameOut = ['melModel_comp_', 'norm-', normStyle{i}, '_', linLog{k},  '.png'];
                            export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)                
                        end
                
                % residual plot
                fig = figure('Color','w');
                
                    plotResiduals(brainardStats, thapanStats, najjarStats, titleStr, handles);
                    
                        % autosave the figure      
                        if style.imgOutautoSavePlot == 1
                            fileNameOut = ['melModel_residuals_', 'norm-', normStyle{i}, '_', linLog{k},  '.png'];
                            export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)                
                        end
            
            end
            
        end