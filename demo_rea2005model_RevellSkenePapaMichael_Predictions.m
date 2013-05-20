% Demonstrates the Rea et al. 2005 model for melatonin suppression
function demo_rea2005model_RevellSkenePapaMichael_Predictions()

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
        criteriaVector = 300; % 300 is the default value
        % FWHM = 10; % of the monochromatic light spectra
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
    
        tic;    
        lambdaSimSPD = [437 479 627 532];
            % 479 and 627 nm from Papamichael et al. (2012), with 10 nm hbw
                % http://dx.doi.org/10.1177/0748730411431447
            % 437, 479 and 532 from Revell et al. (2010)
                % http://dx.doi.org/10.3109/07420528.2010.516048
        FWHM = [10 10 10 10]; % of the monochromatic light spectra

        % create monochromatic light, irradiance weighed
        [SPD.raw{1}, SPD.orig{1}, irrad] = createMonochromaticLightSPDs(lambda, lambdaSimSPD, FWHM, photDensityCnst, irradMin, irradMax, irradRes, path);

        % correct for ocular media transmission
        % [SPD.filt, lensFilter] = correctForOcularMedia(lambdaSimSPD, SPD, irrad, age, lambda, offset, path);             
                    
    %% Define the experimental lights
    
        % From Revell et al. (2010), 6 conditions (Table 1)
            % 1) 437 nm     4.2x10e13   19.1uW/cm2
            % 2) 479 nm     2.5x10e13   10.4uW/cm2
            % 3) 532 nm     6.4x10e13   23.8uW/cm2
            % 4) 437+479 nm     6.7x10e13   29.5uW/cm2
            % 5) 479+479 nm     5.0x10e13   20.8uW/cm2
            % 6) 479+532 nm     8.9x10e13   34.2uW/cm2
            
        % From Papamichael et al. (2012), 7 + 7 conditions (Table 1)
            % 7) 479 nm     1x10e13
            % 8) 479 nm     5x10e13
            % 9) 479 nm     1x10e14
            % 10) 627 nm     5x10e13            
            % 11) 479+627 nm     6x10e13
            % 12) 479+627 nm     1x10e14
            % 13) 479+627 nm     1.5x10e14
            
            % 14) 479 nm     2.5x10e13
            % 15) 627 nm     5x10e13
            % 16) 627 nm     1x10e14
            % 17) 627 nm     3x10e14
            % 18) 479+627 nm     7.5x10e13
            % 19) 479+627 nm     1.3x10e14
            % 20) 479+627 nm     3.3x10e14
            
        % pick the correct wavelength (created monochromatic SPDs), the
        % correct photon density, and the corresponding melatonin
        % suppression (read from a table)
        
            fileName = 'revellSkenePapaMichael_lightTable.txt';
            fid = fopen(fullfile(path.dataIn,fileName));
            protocolSpecs = textscan(fid, '%s %n %n %n %n %n %n %n %s', 'HeaderLines', 1, 'Delimiter', '\t');
        
            for i = 1 : length(protocolSpecs{2})
                wavelengthInd{i}(1) = protocolSpecs{2}(i);
                if ~isnan(protocolSpecs{3}(i))
                    wavelengthInd{i}(2) = protocolSpecs{3}(i);
                end                    
            end
            
            for i= 1 : length(protocolSpecs{4})
                photonDensity{i}(1) = protocolSpecs{4}(i);
                if ~isnan(protocolSpecs{5}(i))
                    photonDensity{i}(2) = protocolSpecs{5}(i);
                end
            end
                
            melSuppressionValues = protocolSpecs{7} / 100;
            melSuppressionSD = protocolSpecs{8};
    
            
            % create the xTic labels
            for i = 1 : length(protocolSpecs{1})
               phDensStr = num2str(sum(photonDensity{i}), '%1.2E')
               xTickLabels{i} = sprintf('%s%s%s%s%s%s%s', ' ', num2str(i), ') ', protocolSpecs{1}{i}, 'nm, ', phDensStr, ' ph/Cm^2/sec');
            end
            
        % NOW Create the final experimental lights with three different
        % ocular media filtering
        for i = 1 : length(wavelengthInd)
            SPD.experim{1}(:,i) = createExperimentalLight(lambdaSimSPD, lambda, SPD.raw, wavelengthInd{i}, photonDensity{i});            
            % plot(lambda, SPD.experim{1}(:,i)); pause(0.6);            
        end
                
        % correct for ocular media transmission
        
            cd(path.ocularMedia)
            lensFilter{1} = agedLensFilter(age(1), lambda, offset); % 25 yr std observer
            lensFilter{2} = agedLensFilter(age(2), lambda, offset); % 65 yr std observer
        
            for i = 1 : length(wavelengthInd)
                SPD.experim{2}(:,i) = SPD.experim{1}(:,i) .* lensFilter{1};
                SPD.experim{3}(:,i) = SPD.experim{1}(:,i) .* lensFilter{2};
            end
            cd(path.mainCode)
            
        
    %% Compute responses and parameters
    
        S_cornea.Vl(isnan(S_cornea.Vl)) = 0;

        cd(path.photoreception)        
        for i = 1 : length(SPD.experim)
            for j = 1 : length(wavelengthInd)                
                
                CLA{i}(j) = CLAfuncComp(SPD.experim{i}(:,j), S_cornea, S_retina, spectralCrossover, modelType, a_cone);
                CS{i}(j)  = .7 * (1 - (1./(1 + (CLA{i}(j)/355.7) .^ (1.1026))));
                
                illum{i}(j) = 683 * trapz(S_cornea.Vl .* SPD.experim{i}(:,j)); % lux
                irradiance{i}(j) = sum(SPD.experim{i}(:,j)) * 10^2; % uW/cm2
            end
        end
        
        cd(path.mainCode)
        timing = toc
            
    %% PLOT The Results
    
        x = 1 : 1 : length(wavelengthInd);
    
        fig = figure('Color', 'w')
            set(fig, 'Position', [0.01*scrsz(3) 0.21*scrsz(4) 0.78*scrsz(3) 0.7*scrsz(4)])
        
            rows = 3;
            cols = 3;
        
            i = 1;
            sp(i) = subplot(rows,cols,i);
            p1 = plot(x, CS{1}, 'or', x, CS{2}, 'og', x, CS{3}, 'ob');
                lab(i,1) = ylabel('CS'); lab(i,2) = xlabel('Light Condition');
            
            i = 2;
            sp(i) = subplot(rows,cols,4);
            p2 = semilogy(x, illum{1}, 'or', x, illum{2}, 'og', x, illum{3}, 'ob');
                lab(i,1) = ylabel('Illuminance [lux]'); lab(i,2) = xlabel('Light Condition');
                legend('Cornea', '25 yr lens', '65 yr lens', 3, 'Location', 'NorthWest')
                legend('boxoff')
            
            i = 3;
            sp(i) = subplot(rows,cols,7);
            p3 = semilogy(x, irradiance{1}, 'or', x, irradiance{2}, 'og', x, irradiance{3}, 'ob');
                lab(i,1) = ylabel('Irradiance [uW/cm^2]'); lab(i,2) = xlabel('Light Condition');
            
            i = 4;
            sp(i) = subplot(rows,cols,[2 3 5 6]);
            p4 = plot(x, CS{1}, 'og', x, melSuppressionValues, 'ok', x, abs(CS{1}' - melSuppressionValues), 'oy');
                lab(i,1) = ylabel('CS, Melatonin Suppression'); lab(i,2) = xlabel('  ');
                
                set(gca, 'XTick', x, 'XTickLabel', '')
                t = text(x, zeros(length(x),1), xTickLabels, 'Rotation', -90);
                
                legend('CS (cornea)', 'Mel.Suppression', 'Residuals (abs.)', 3, 'Location', 'NorthWest')
                legend('boxoff')
            
                
        % STYLE
            markerBaseSize = 4;
            set([p1 p2 p3], 'MarkerSize', markerBaseSize)
            set([p4], 'MarkerSize', markerBaseSize+3)
            set(p4(1), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
            set(p4(2), 'MarkerFaceColor', 'k')
            set(p4(3), 'MarkerFaceColor', [0.871 0.490 0], 'MarkerEdgeColor', [0.871 0.490 0])
            % set([p1 p2 p3 p4], 'MarkerFaceColor', 'b', 'MarkerSize', markerBaseSize+2)
            
            set(sp, 'XLim', [0.5 length(wavelengthInd)+0.5], 'XTick', x)
            
            set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)
            set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')
            set(t, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')
                
        % autosave the figure [PNG]       
        fileNameOutBase = 'melPred_revellSkenePapaMichael_reaModel';
        
        if style.imgOutautoSavePlot == 1
            fileNameOut = [fileNameOutBase, '.png'];
            try
                export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)                
            catch err
                err
                warning('%s\n%s\n%s\n%s', err.identifier, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"', ...
                      '"Home" -> "Environment" -> "Set Path" in Matlab2012b and newer with the updated GUI')   
            end
        end     
        
        % autosave the figure [SVG]        
        if style.imgOutautoSavePlot == 1
            fileNameOut = [fileNameOutBase, '.svg'];
            try
                plot2svg((fullfile(path.figuresOut, fileNameOut)))                
            catch err
                warning('%s\n%s\n%s\n%s', err.identifier, 'Figure not saved most likely because you have not installed plot2svg from Matlab File Exchange!', ...
                      'Download it from: https://www.mathworks.com/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures, and "File -> Set Path -> Add Folder"', ...
                      '"Home" -> "Environment" -> "Set Path" in Matlab2012b and newer with the updated GUI')   
            end
        end       