function [brainardStats, thapanStats, najjarStats] = plotModelComparison(j, fig, scrsz, style, x, y, normStyle, linLog, yStr, titleStr, path, handles)

    %% Import Melatonin data

        % Brainard and Thapan Data
        % [wavelength  response BrainardYES]
        BandT = [420 0.256 1; ...
                 424 0.815 0; ...
                 440 0.953 1; ...
                 456 1.000 0; ...
                 460 1.000 1; ...
                 472 0.8560 0; ...
                 480 0.7259 1; ...
                 494 0.5869 0; ...
                 505 0.7916 1; ...
                 520 0.5202 0; ...
                 530 0.3958 1; ...
                 548 0.1428 0; ...
                 555 0.1089 1; ...
                 575 0.0554 1; ...
                 600 0.0282 1];

                waveBT = BandT(:,1);
                BT = BandT(:,2);

                booleanBrainard = logical(BandT(:,3));
                waveBrainard = waveBT(booleanBrainard);
                Brainard =  BT(booleanBrainard);
                waveThapan =  waveBT(~booleanBrainard);
                Thapan =  BT(~booleanBrainard);

        % NAJJAR et al.
        najjar = [420.0000    0.1763    0.1346; ...
                  440.0000    0.7338    0.1588; ...
                  460.0000    0.9796    0.1794; ...
                  480.0000    0.9406    0.1771; ...
                  500.0000    1.0000    0.1252; ...
                  530.0000    0.9725    0.1422; ...
                  560.0000    0.4260    0.1767; ...
                  590.0000    0.0743    0.2069; ...
                  620.0000    0.1676    0.1719];

    
    %% Change figure size
    
        if strcmp(linLog, 'lin') == 1
            set(fig, 'Position', [0.01*scrsz(3) 0.05*scrsz(4) 0.91*scrsz(3) 0.7*scrsz(4)])
        elseif strcmp(linLog, 'log') == 1
            set(fig, 'Position', [0.03*scrsz(3) 0.02*scrsz(4) 0.91*scrsz(3) 0.7*scrsz(4)])
        else
           error('linLog definition incorrect'); 
        end

        % Set correct legend labels
        
            ocuMediaStrLabels = {'Aphakic', '25yr.', '65yr.'};
    
    %% DEFINE HOW TO NORMALIZE THE DATA
    
        disp(['Normalization style: ', normStyle])

        % get normalization factors
        if strcmp(normStyle, 'toUnity')

            offset = 40; % ignore short-wavelength horrors that arise for aphakic ones
            norm{1} = max(y{1}(offset:end));
            norm{2} = max(y{2}(offset:end));
            norm{3} = max(y{3}(offset:end));

            % Scaler, for Najjar et al.
            scaler = 1;        

        elseif strcmp(normStyle, 'to25years')

            offset = 40; % ignore short-wavelength horrors that arise for aphakic ones
            if strcmp(yStr, 'Rea/LRC')
                norm{1} = 1 / 3.0653;
                norm{2} = 1;
                norm{3} = 1 / 0.6161;
                % taken manually from demo_rea2005model_monochromaticSpectra
            else
                norm{1} = max(y{2}(offset:end));
                norm{2} = max(y{2}(offset:end));
                norm{3} = max(y{2}(offset:end));
            end

            % Scaler, for Najjar et al.
            atWavelength = 460;
            wavelIndex = find(x == atWavelength);

            if ~isempty(wavelIndex)                   
                normY = y{3}/norm{3};
                scaler = normY(wavelIndex);
            else
                warning('Unscaled Najjar points, too coarse wavelength definition?')
                scaler = 1;
            end        

        elseif strcmp(normStyle, 'toLongWavelength')

            offset = 40; % ignore short-wavelength horrors that arise for aphakic ones
            if strcmp(yStr, 'Rea/LRC')
                norm{1} = 1 / 3.0653;
                norm{2} = 1;
                norm{3} = 1 / 0.6161;
                % taken manually from demo_rea2005model_monochromaticSpectra
            else
                norm{1} = max(y{2}(offset:end));
                norm{2} = max(y{2}(offset:end));
                norm{3} = max(y{2}(offset:end));
            end

            % Scaler, for Najjar et al.
            atWavelength = 530;
            wavelIndex = find(x == atWavelength);

            if ~isempty(wavelIndex)      
                normY = y{3}/norm{3};
                scaler = normY(wavelIndex);
            else
                warning('Unscaled Najjar points, too coarse wavelength definition?')
                scaler = 1;
            end        

        else

            errordlg('How did you want to normalize the data? typo in normStyle?')

        end
    
    %% CALCULATE STATS
    
        % get the data corresponding to the different data sets
        
        for i = 1 : length(waveBrainard)
            try
                brainardIndices(i,1) = find(waveBrainard(i) == x);
            catch err
                warning('Too coarse spectral resolution, no match found')
                brainardIndices(i,1) = NaN;
            end
        end
        for i = 1 : length(waveThapan)
            try
                thapanIndices(i,1)   = find(waveThapan(i) == x);
            catch err
                warning('Too coarse spectral resolution, no match found')
                thapanIndices(i,1) = NaN;
            end
        end
        for i = 1 : length(najjar(:,1))
            try
                najjarIndices(i,1)   = find(najjar(i,1) == x);
            catch err
                warning('Too coarse spectral resolution, no match found')
                najjarIndices(i,1) = NaN;
            end
        end
        
        % use a subfunction to calculate the fitting stats        
            
            % 25 yr old Std Observer
            yFit_25yrs   = (y{2} / norm{2});            
                [rows,cols] = size(yFit_25yrs); % transpose check
                if rows < cols; yFit_25yrs = yFit_25yrs'; end      

            % 65 yr old Std Observer
            yFit_elderly = (y{3} / norm{3});            
                [rows,cols] = size(yFit_elderly); % transpose check
                if rows < cols; yFit_elderly = yFit_elderly'; end
            
            % Debug
            wavel = x(najjarIndices); 
            najjarFit = yFit_elderly(najjarIndices);   
            if j == 1
                scaler
                [wavel najjarFit najjar(:,1) (scaler*najjar(:,2))]
            end
            
            % calculate stats
            if strcmp(linLog, 'lin') == 1
                
                brainardStats.R2 = rsquare(Brainard, yFit_25yrs(brainardIndices));
                thapanStats.R2 = rsquare(Thapan, yFit_25yrs(thapanIndices));
                najjarStats.R2 = rsquare(scaler*najjar(:,2), yFit_elderly(najjarIndices));
                
            elseif strcmp(linLog, 'log') == 1
                
                brainardStats.R2 = rsquare(log10(Brainard), log10(yFit_25yrs(brainardIndices)));
                thapanStats.R2 = rsquare(log10(Thapan), log10(yFit_25yrs(thapanIndices)));
                najjarStats.R2 = rsquare(log10(scaler*najjar(:,2)), log10(yFit_elderly(najjarIndices)));
                
            end
            
            % CALCULATE RESIDUALS    
            brainardStats.residuals.y = Brainard - yFit_25yrs(brainardIndices);
            brainardStats.residuals.x = waveBrainard;
            thapanStats.residuals.y   = Thapan - yFit_25yrs(thapanIndices);
            thapanStats.residuals.x   = waveThapan;
            najjarStats.residuals.y   = (scaler*najjar(:,2)) - yFit_elderly(najjarIndices);
            najjarStats.residuals.x   = najjar(:,1);            
               
    %% CALCULATE PEAK SHIFTS
        
        offset = 20;
    
        [val,ind] = max(y{1}(offset:end) / norm{1});
            peakNm(1) = x(offset+ind);
        [val,ind] = max(y{2}(offset:end) / norm{2});
            peakNm(2) = x(offset+ind);
        [val,ind] = max(y{3}(offset:end) / norm{3});
            peakNm(3) = x(offset+ind);           
    
        
    %% PLOT
    
        % remove the second column if there is any
        y{1} = y{1}(:,1);
        y{2} = y{2}(:,1);
        y{3} = y{3}(:,1);
       
        % The model predictions (Aphakic / 25yr / 65 yr)
        if strcmp(linLog, 'lin') == 1
            p1 = plot(x, y{1} / norm{1}, ...
                             x, y{2} / norm{2}, ...
                             x, y{3} / norm{3})
            
        elseif strcmp(linLog, 'log') == 1
            p1 = plot(x, log10(y{1} / norm{1}), ...
                            x, log10(y{2} / norm{2}), ...
                            x, log10(y{3} / norm{3}))
                      
        else
            error('linLog definition incorrect');
        end
        
            hold on

             % Brainard & Thapan data         

                if strcmp(linLog, 'lin') == 1
                    p1b(1) = plot(waveBrainard, Brainard, 'bo');
                    p1b(2) = plot(waveThapan, Thapan, 'ro');
                elseif strcmp(linLog, 'log') == 1
                    p1b(1) = plot(waveBrainard, log10(Brainard), 'bo');
                    p1b(2) = plot(waveThapan, log10(Thapan), 'ro');
                else
                   error('linLog definition incorrect'); 
                end
                set(p1b(1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
                set(p1b(2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')

            % Najjar

                if strcmp(linLog, 'lin') == 1
                    p1b(3) = errorbar(najjar(:,1), scaler*najjar(:,2), najjar(:,3), 'ko');
                elseif strcmp(linLog, 'log') == 1
                    logError_LO = log10(scaler*najjar(:,2) - najjar(:,3));
                    logError_HI = log10(scaler*najjar(:,2) + najjar(:,3));
                    p1b(3) = errorbar(najjar(:,1), log10(scaler*najjar(:,2)), logError_LO, logError_HI, 'ko');
                else
                    error('linLog definition incorrect');
                end
                set(p1b(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
                hold off 

    %% Style
    
        r2prec = '%1.2f';
        peakPrec = '%3.0f';
         
        leg(1) = legend(['Aphakic, \lambda_m_a_x = ', num2str(peakNm(1), peakPrec), ' nm'],...
                        ['25 yrs., \lambda_m_a_x = ', num2str(peakNm(2), peakPrec), ' nm'],...
                        ['65 yrs., \lambda_m_a_x = ', num2str(peakNm(3), peakPrec), ' nm'],...            
                        ['Brainard, R^2 = ', num2str(brainardStats.R2,r2prec)],...
                        ['Thapan, R^2 = ', num2str(thapanStats.R2,r2prec)],...
                        ['Najjar, R^2 = ', num2str(najjarStats.R2,r2prec)]);
        legend('boxoff')
        
        lab(1) = xlabel('Wavelength [nm]');     
        lab(2) = ylabel(yStr);
        tit = title(titleStr); 
        
        % line plots            
        set(p1(1), 'Color', style.colorPlot(1,:))
        set(p1(2), 'Color', style.colorPlot(2,:))
        set(p1(3), 'Color', style.colorPlot(3,:))

        % legend % text
        set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)        
        set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')    
        set(tit, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')
        
        % subplots (gca)
        set(gca, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)   
        set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)        
        set(gca, 'XLim', [400 650])
        if strcmp(linLog, 'lin') == 1
            set(gca, 'YLim', [0 1.2])
            if strcmp(normStyle, 'to25years') || strcmp(normStyle, 'to25years')
                set(gca, 'YLim', [0 2])
            end
        elseif strcmp(linLog, 'log') == 1
            set(gca, 'YLim', [-2.5 0.6])
            set(leg, 'Location', 'South')
        end
