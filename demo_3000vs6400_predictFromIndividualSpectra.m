function demo_3000vs6400_predictFromIndividualSpectra()

    % Predicts the melatonin suppression for the goggle study

    % Petteri Teikari, 2013, LRC, Troy, NY, USA
    % teikap@rpi.edu
    close all    
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    
    %% SETTINGS

        % FOLDERS
        mFilename = 'demo_3000vs6400_predictFromIndividualSpectra'; 
        pathF = mfilename('fullpath');
        path.mainCode = strrep(pathF, mFilename, '');

        % FOLDERS
        path = initFolders();  
            
         % Styling
        style.fontName = 'latin modern roman';
        style.fontBaseSize = 9;
        style.fontLabelWeight = 'bold';
       
        % FIGURE autosave options
        style.imgOutRes       = '-r300'; % dpi
        style.imgOutAntiAlias = '-a1';   % '-a0' least, '-a4' maximum
        style.imgOutautoSavePlot = 1;
        
        % MODEL Parameters
        spectralCrossover = 'sharp'; 
            % spectralCrossover = 'smooth'; 
        modelType = 'original';  
            % modelType = 'origWithCones';
        a_cone = 0.0;
        
    %% IMPORT THE Spectral Sensitivities of Photoreceptors        
    
        % Corneal Sensitivities
        cd(path.photoreception)
        S_cornea = import_CornealSensitivities(path.templates);
        cd(path.mainCode)
        
        % Retinal Sensitivities
        cd(path.photoreception)
        S_retina = import_RetinalSensitivities(path.nomogram);
        cd(path.mainCode)
        
    %% IMPORT MELATONIN SUPPRESSION DATA
    
        % define the filenames for MELATONIN SUPPRESSION
        fileNames = {'MelatoninSuppression3000K_1.csv'; 'MelatoninSuppression3000K_2.csv'; 'MelatoninSuppression6400K.csv'};        
            headerLines = 3;
            delimiterIn = ',';
            dataIn = cell(length(fileNames),1);
        
            % import the data
            for i = 1 : length(fileNames)            

               dataRaw = importdata(fullfile(path.dataIn, fileNames{i}), delimiterIn, headerLines);

               % split the structure to variables
               dataIn{i}.melSuppression = dataRaw.data(:,2:end);
               dataIn{i}.subjects = dataRaw.data(:,1);
               dataIn{i}.colheaders = dataRaw.colheaders;
               dataIn{i}.illuminances = dataIn{i}.colheaders; 

                    dataIn{i}.illuminances = strrep(dataIn{i}.illuminances, 'sub #', '');
                    dataIn{i}.illuminances = strrep(dataIn{i}.illuminances, ' lux', '');
                    dataIn{i}.illuminances = str2double(dataIn{i}.illuminances);
                    dataIn{i}.illuminances = dataIn{i}.illuminances(2:end); % get rid of subject column

               dataIn{i}.description = cell2mat(dataRaw.textdata(1,1));
                    dataIn{i}.description = strrep(dataIn{i}.description, ',', '');

               CCT_raw = dataRaw.textdata(2,1);

                    fields = textscan(cell2mat(CCT_raw), '%s%s%s%s', 'Delimiter', ',');
                    dataIn{i}.CCT = str2double(strrep(fields{2}, 'K', ''));
            end
         
    %% IMPORT the light SPD
            
        % define the filenames for the LIGHT SPD
        fileNamesSPD = {'Average spd for 3000K experiment.txt' ; '6400K_spd.txt'};
            headerLines = 1;
            delimiterIn = '\t';
            lambdaNEW = (380:1:780)';
            
            % import the tape that was in front of the light sources that
            % filtered the SPD slightly
            fileNamesFilter = 'ongoggle_tapetrans.txt';
                raw = importdata(fullfile(path.dataIn, fileNamesFilter), delimiterIn, headerLines);
                tapeFilter = interp1(raw.data(:,1), raw.data(:,2), lambdaNEW, 'linear');
                    tapeFilter(isnan(tapeFilter)) = 1; % pad with ONEs 
            
            
            % import the data
            for i = 1 : length(fileNamesSPD)                
                dataRaw = importdata(fullfile(path.dataIn, fileNamesSPD{i}), delimiterIn, headerLines);
                
                wavelength = dataRaw.data(:,1);
                irradiance = dataRaw.data(:,2);
                
                irradianceNEW = interp1(wavelength, irradiance, lambdaNEW, 'linear');
                    %     sp(i) = subplot(1,2,i);
                    %         plot(wavelength, irradiance, 'r', lambdaNEW, irradianceNEW, 'k')                
                irradiance = irradianceNEW;
                    irradiance(isnan(irradiance)) = 0; % get rid of NaNs
                    irradiance = irradiance .* tapeFilter; % correct for the tape filter
                wavelength = lambdaNEW;             
                
                if i == 1 % 3000 K
                    
                    ind = 1;
                    dataIn{ind}.wavelength = wavelength;
                    dataIn{ind}.irradiance = irradiance;                    
                    
                    ind = 2;
                    dataIn{ind}.wavelength = wavelength;
                    dataIn{ind}.irradiance = irradiance;                
                    
                elseif i == 2 % 6500/6400 K
                    
                    ind = 3;
                    dataIn{ind}.wavelength = wavelength;
                    dataIn{ind}.irradiance = irradiance;
                    
                end
                
            end       
            
    %% CORRECT THE LIGHT SPECTRUM to have to correct ILLUMINANCE
    
        % import V(lambda)
        dataRaw = importdata(fullfile(path.dataIn, 'v_lambda_linCIE2008v10e_LIN_380to780nm_1nm.txt'), '\t', 1);
            vLambda = dataRaw.data(:,2);
            vLambda(isnan(vLambda)) = 0;
            
        disp(' ');
        for i = 1 : length(dataIn)
            for j = 1 : length(dataIn{i}.illuminances)
                
                irradIN = dataIn{i}.irradiance;
                illuminanceTARGET = dataIn{i}.illuminances(j);
                
                luxCalc = 683 * trapz(vLambda .* irradIN);
                ratioLux = illuminanceTARGET / luxCalc;
                irradOUT = irradIN * ratioLux
                
                luxCalcCheckUp =  683 * trapz(vLambda .* irradOUT);
                disp(['i=', num2str(i), ', j=', num2str(j), ', target=', num2str(illuminanceTARGET), ', luxCalc=', num2str(luxCalc), ', ratio=', num2str(ratioLux), ', check=', num2str(luxCalcCheckUp)])
                                
                dataIn{i}.irradianceMatrix_woPupil(:,j) = irradOUT;
                
            end
        end
        
    %% IMPORT THE PUPIL DIAMETERS
    
        % define the filenames for the LIGHT SPD
        fileNamesPupil = {'pupilDiameters_3000K_1.txt'; 'pupilDiameters_3000K_2.txt'; 'pupilDiameters_6400K.txt'};
            headerLines = 3;
            delimiterIn = '\t';
            lambdaNEW = (380:1:780)';
            
            for i = 1 : length(fileNamesPupil)
                raw = importdata(fullfile(path.dataIn, fileNamesPupil{i}), delimiterIn, headerLines);
                for j = 1 : length(dataIn{i}.illuminances)
                    if i == 1 % with the DARK
                        dataIn{i}.pupilDiameter(j) = raw.data(j+1,2);
                    else
                        dataIn{i}.pupilDiameter(j) = raw.data(j,2);
                    end
                end
                % a = dataIn{i}.pupilDiameter;
            end
            
    %% Correct the illuminances with the pupil diameters in relation to the biggest pupil diameter
    
        % the biggest pupil diameter
        referenceDiameter = max(dataIn{1}.pupilDiameter);
        referenceArea = (referenceDiameter/2)^2 * pi;        
        
        for i = 1 : length(dataIn)
            for j = 1 : length(dataIn{i}.illuminances)                
                pupilArea = (dataIn{i}.pupilDiameter(j) / 2)^2 * pi;
                dataIn{i}.pupilAreaRelative(j) = pupilArea / referenceArea;
                
                dataIn{i}.irradianceMatrix(:,j) = dataIn{i}.irradianceMatrix_woPupil(:,j) * dataIn{i}.pupilAreaRelative(j);
                
                dataIn{i}.retinalIlluminance(j) = 683 * trapz(vLambda .* dataIn{i}.irradianceMatrix(:,j));                
            end
        end
            
    %% CALCULATE THE CS VALUES
        
        cd(path.photoreception)
        for i = 1 : length(dataIn)
            for j = 1 : length(dataIn{i}.illuminances)                 
                [dataIn{i}.CLA(j), dataIn{i}.Mel(j), dataIn{i}.SWS(j), dataIn{i}.Vl(j), dataIn{i}.Cones(j), dataIn{i}.Rod(j), dataIn{i}.CLAComp(j), dataIn{i}.oppFlag(j)] = ...
                    CLAfuncComp(dataIn{i}.irradianceMatrix(:,j), S_cornea, S_retina, spectralCrossover, modelType, a_cone);
                dataIn{i}.CS(j)  = .7 * (1 - (1./(1 + (dataIn{i}.CLA(j)/355.7).^(1.1026)))); % CSCalc_postBerlin_12Aug2011
            end
            
            % aa = dataIn{i}.Rod % are NaN with melanopsin-only response
        end
        cd(path.mainCode)
        
        
    %% PLOT    
    
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.01*scrsz(3) 0.21*scrsz(4) 0.78*scrsz(3) 0.7*scrsz(4)])
        
        % subplot layout
        rows = 2;
        cols = length(fileNames);
        
        % MELATONIN SUPPRESSION
    
            for i = 1 : length(fileNames) 
                sp(i) = subplot(rows,cols,i);
                    [dataIn{i}.melSupprMean, dataIn{i}.melSupprSD] = plotSuppressionValues(fileNames{i}, dataIn{i}.subjects, dataIn{i}.melSuppression, dataIn{i}.colheaders, dataIn{i}.illuminances, dataIn{i}.CCT, dataIn{i}.description, style);
            end

        % LIGHT SPD
            i = 4;
            sp(i) = subplot(rows,cols,i);
            
                pSPD = plot(wavelength, dataIn{1}.irradiance, 'r', wavelength, dataIn{3}.irradiance/max(dataIn{3}.irradiance), 'b', wavelength, vLambda, 'g');
                    leg(1) = legend('3000 K', '6400 K', 'V(\lambda)');
                        legend('boxoff')
                        xlim([380 780])
                        ylim([0 1.05])
        
        % SIGMOID (3000K)
            i = 5;
            sp(i) = subplot(rows,cols,[i i+1]);            
            
                hold on 
                pSigmoid1(1) = errorbar(dataIn{1}.CLA, dataIn{1}.melSupprMean, dataIn{1}.melSupprSD, 'ko');
                pSigmoid1(2) = errorbar(dataIn{2}.CLA, dataIn{2}.melSupprMean, dataIn{2}.melSupprSD, 'ko');
                    set(get(pSigmoid1(1),'Parent'), 'XScale', 'log')                   
                
                pSigmoid2(1) = errorbar(dataIn{3}.CLA, dataIn{3}.melSupprMean, dataIn{3}.melSupprSD, 'ko');
                    set(get(pSigmoid2(1),'Parent'), 'XScale', 'log')

                    % plot the CS (sigmoid-like curve)
                    minX = 1; maxX = 10000;
                    CS_x = (minX : 1 : maxX)';
                    CS_vector = zeros(length(CS_x),1);
                    for i = minX : maxX
                       CS_vector(i) = .7 * (1 - (1./(1 + (i/355.7).^(1.1026)))); % CSCalc_postBerlin_12Aug2011
                    end
                    
                    pSigmoidCS = plot(CS_x, CS_vector, 'k');
                    
                
                    
                    
                % style
                leg(2) = legend('3000 K_a', '3000 K_b', '6400 K', 'CS', 4, 'Location', 'NorthWest');
                    legend('boxoff')
                
                xlim([1 10000])
                ylim([-0.3 0.7])
                
                lab(1) = ylabel('Melatonin Suppression');
                lab(2) = xlabel('CL_a');
                
                markerBaseSize = 6;
                set(pSigmoid1, 'MarkerFaceColor', 'r', 'MarkerSize', markerBaseSize, 'Color', [.3 .3 .3])
                set(pSigmoid1(2),'MarkerFaceColor',[1 0.60 0.78])
                set(pSigmoid2, 'MarkerFaceColor', 'b', 'MarkerSize', markerBaseSize, 'Color', [.3 .3 .3])
                
                disp(' ');
                disp('The Subplot 5 data:')
                disp(' ... '); disp(' ');
                
                disp('3000 K'); disp(' ');
                dataMatrix_headers = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', 'Illuminance', 'CLA', 'Melatonin Suppr. Mean', 'Melatonin Suppr. SD', 'CS', 'Pupil Diameter [mm]', 'Relative Pupil Area', 'Retinal Illuminance');
                    disp(dataMatrix_headers)
                    
                    i = 1;
                    for j = 1 : length(dataIn{i}.illuminances)                 
                        dataMatrix{1}(j,:) = [dataIn{i}.illuminances(j) dataIn{i}.CLA(j) dataIn{i}.melSupprMean(j) dataIn{i}.melSupprSD(j) dataIn{i}.CS(j) dataIn{i}.pupilDiameter(j) dataIn{i}.pupilAreaRelative(j) dataIn{i}.retinalIlluminance(j)];                        
                    end
                    
                    i = 2;
                    for jk = 1 : length(dataIn{i}.illuminances)                 
                        dataMatrix{1}(j+jk,:) = [dataIn{i}.illuminances(jk) dataIn{i}.CLA(jk) dataIn{i}.melSupprMean(jk) dataIn{i}.melSupprSD(jk) dataIn{i}.CS(jk) dataIn{i}.pupilDiameter(jk) dataIn{i}.pupilAreaRelative(jk) dataIn{i}.retinalIlluminance(jk)];
                    end
                    disp(dataMatrix{1})
                    
                disp(' '); disp('6400 K'); disp(' ');                                    
                    
                    i = 3;
                    for j = 1 : length(dataIn{i}.illuminances)                 
                        dataMatrix{2}(j,:) = [dataIn{i}.illuminances(j) dataIn{i}.CLA(j) dataIn{i}.melSupprMean(j) dataIn{i}.melSupprSD(j) dataIn{i}.CS(j) dataIn{i}.pupilDiameter(j)  dataIn{i}.pupilAreaRelative(j) dataIn{i}.retinalIlluminance(j)];
                    end
                    disp(dataMatrix{2})
                
                % write to disk
                fileOut = 'melPred_PolyChrom_3000K_vs_6400K_withTape_w_PupilMeasAsIlluminanceCorrection.txt';
                fileWithPath = fullfile(path.dataOut, fileOut);                
                dlmwrite(fileWithPath, dataMatrix_headers, 'delimiter', '')
                dlmwrite(fileWithPath, '3000K', '-append', 'delimiter', '')
                dlmwrite(fileWithPath, dataMatrix{1}, '-append', 'delimiter', '\t')
                dlmwrite(fileWithPath, '6400K', '-append', 'delimiter', '')
                dlmwrite(fileWithPath, dataMatrix{2}, '-append', 'delimiter', '\t')
                    
        % Something
        %             i = 6;
        %             sp(i) = subplot(rows,cols,i);      
        set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')
            
        % autosave the figure [PNG]        
        if style.imgOutautoSavePlot == 1
            fileNameOut = ['melPred_PolyChrom_3000K_vs_6400K_withTape_w_PupilMeasAsIlluminanceCorrection', '.png'];
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
            fileNameOut = ['melPred_PolyChrom_3000K_vs_6400K_withTape_w_PupilMeasAsIlluminanceCorrection', '.svg'];
            try
                plot2svg((fullfile(path.figuresOut, fileNameOut)))                
            catch err
                err
                warning('%s\n%s\n%s\n%s', err.identifier, 'Figure not saved most likely because you have not installed plot2svg from Matlab File Exchange!', ...
                      'Download it from: https://www.mathworks.com/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures, and "File -> Set Path -> Add Folder"', ...
                      '"Home" -> "Environment" -> "Set Path" in Matlab2012b and newer with the updated GUI')   
            end
        end     
    
        
    %% SUBFUNCTIONS
    
        function [meanValue, sdValue] = plotSuppressionValues(fileName, subjects, melSuppression, colheaders, illuminances, CCT, description, style)
            
            % calculate means and SDs
            meanValue = nanmean(melSuppression);
            sdValue   = nanstd(melSuppression);
            
            hold on
            
                offset = 0.08;
                x_init = [1 2 3];
            
                % individual suppression values
                x_indiv = x_init - offset;
                p1 = plot(x_indiv, melSuppression', '*');
                        
                % the mean of condition                
                x_mean = x_init + offset;
                p2 = errorbar(x_mean, meanValue, sdValue, 'ko');
                
                % annotate the mean value to the plot
                annotOffset = 0.15;
                prec = '%1.2f';
                for ij = 1 : length(meanValue)
                    textValues{ij} = sprintf('%s%s%s', num2str(meanValue(ij),prec), '\pm', num2str(sdValue(ij),prec));
                    t(ij) = text(x_init(ij)+annotOffset, meanValue(ij), textValues{ij});
                end                
                
                
            hold off
            
            % style
            
                luxStr = repmat(' lux', length(illuminances),1);            
                set(gca, 'XTick', x_init, 'XTickLabel', [num2str(illuminances') luxStr])

                titleStr = sprintf('%s\n%s', description, [num2str(CCT), ' K']);
                tit = title(titleStr);

                lab = ylabel('Melatonin Suppression');
                
                set(gca, 'XLim', [x_indiv(1)-(2*offset) x_indiv(end)+(2*offset)], 'YLim', [-1 1])
                
                set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')  
                set(tit, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')
                set(t, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)
                
                markerBaseSize = 5;
                set(p1, 'MarkerSize', markerBaseSize)
                set(p2, 'MarkerFaceColor', 'b', 'MarkerSize', markerBaseSize+2)
                
