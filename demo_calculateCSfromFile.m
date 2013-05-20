% Demonstrates the Rea et al. 2005 model for melatonin suppression
function demo_calculateCSfromFile()

    % Petteri Teikari, 2013, LRC, Troy, NY, USA
    % teikap@rpi.edu
    close all    
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    
    %% SETTINGS
    
    
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
        
        % Plot switches
        plotON.spectra = 0;
        plotON.components = 1; 
        plotON.criterionComparison = 0;
        
    %% MODEL PARAMETERS
    
        wavelength = (380:1:780)';        
    
        % REA et al. PARAMETERS        
        spectralCrossover = 'sharp'; 
        % spectralCrossover = 'smooth'; 
        modelType = 'original';  
        % modelType = 'origWithCones';
        a_cone = 0.0;
        
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
          
        
    %% DEFINE the MONOCHROMATIC SPECTRA    
    
        % get the filename from the data folder
        [filenameSPD, path.data] = uigetfile('*.csv','Select the Spectrum file (has to be 2 column file with 1 header row, comma-delimited, and a wavelength range between 380-780nm!)');
        
        % import the chosen file     
        tmp = importdata(fullfile(path.data, filenameSPD), ',', 1);        
        
        % assign to different variable names
        headers = tmp.colheaders;
        wavelength_in = tmp.data(:,1); % wavelengths
        SPD_in = tmp.data(:,2:end); % spectral irradiance       
        
        [rows,numberOfSPD] = size(SPD_in);
        
            % quick and dirty fix if the wavelength vector is not 401
            % sample long
            %             wavelength_in
            %             SPD_in, 
            %             wavelength
            for i = 1 : numberOfSPD
                SPD(:,i) = interp1(wavelength_in, SPD_in(:,i), wavelength);
                areNaN = isnan(SPD(:,i));
                SPD(areNaN,i) = 0; % remove NaNs
            end
            
                    
    %% Compute responses and parameters
    
        vLambda = S_cornea. Vl;
        vLambda(isnan(vLambda)) = 0; % NaN -> 0
    
        tic
        cd(path.photoreception)    
        for i = 1 : numberOfSPD
            [CLA(i), Mel, SWS, Vl, Cones, Rod, CLAComp, oppFlag] = CLAfuncComp(SPD(:,i), S_cornea, S_retina, spectralCrossover, modelType, a_cone);
            CS(i)  = .7 * (1 - (1./(1 + (CLA(i)/355.7).^(1.1026)))); % CSCalc_postBerlin_12Aug2011
            illum(i) = 683 * trapz(vLambda .* SPD(:,i));
        end
         
        reaTime = toc;
        
            % rerun with scaled illuminances
            targetIllum = [1 10 500 1000 2500 10000];
            
            for jj = 1 : length(targetIllum)
                
                for i = 1 : numberOfSPD
                    illum = 683 * trapz(vLambda .* SPD(:,i));
                    ratio(jj,i) = targetIllum(jj) / illum;
                    SPD_scaled(:,jj) = ratio(jj,i) * SPD(:,i);
                    [CLAScaled(jj,i), Mel, SWS, Vl, Cones, Rod, CLAComp, oppFlag] = CLAfuncComp(SPD_scaled(:,jj), S_cornea, S_retina, spectralCrossover, modelType, a_cone);
                    CSScaled(jj,i)  = .7 * (1 - (1./(1 + (CLAScaled(jj,i)/355.7).^(1.1026)))); % CSCalc_postBerlin_12Aug2011
                    illumScaled(jj,i) = 683 * trapz(vLambda .* SPD_scaled(:,jj));
                end

            end
            cd(path.mainCode)             
            [ratio(:,1)' CSScaled(:,1)' CS(:,1)']            
            
        % write to disk
        fileOut = 'LightTable_CS_CLA.txt';
        fileWithPath = fullfile(fileOut);    

        dlmwrite(fileWithPath, headers(2:end), 'delimiter', '\t')
        dlmwrite(fileWithPath, CLA, '-append', 'delimiter', '\t')
        dlmwrite(fileWithPath, CS, '-append', 'delimiter', '\t')
        dlmwrite(fileWithPath, illum, '-append', 'delimiter', '\t')
        
            % the scaled
            fileOut = 'LightTable_CS_CLA_scaled.txt';
            fileWithPath = fullfile(fileOut);    

            dlmwrite(fileWithPath, targetIllum, 'delimiter', '\t')
            dlmwrite(fileWithPath, CSScaled(:,1)', '-append', 'delimiter', '\t')
%             

    %% PLOT The Results
    
        figure('Color', 'w')
            p1 = plot(wavelength,SPD); %,'Color',[0.04 0.52 0.78]);
            
            lab(1) = xlabel('Wavelength [nm]');     
            lab(2) = ylabel('W/m^2 ?');    
            xlim([380 780])
           
            set(gca, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)        
            set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')                
            
            leg = legend(headers(2:end));
                legend('boxoff')

            % autosave the figure      
            if style.imgOutautoSavePlot == 1
                fileNameOut = ['lightTableSPD', '.png'];
                try
                    export_fig(fullfile(path.figuresOut, fileNameOut), style.imgOutRes, style.imgOutAntiAlias)
                catch err
                    err.identifier
                    warning('%s\n%s\n%s', err, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"')                        
                end
            end

