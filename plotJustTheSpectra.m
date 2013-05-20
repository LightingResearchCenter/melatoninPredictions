function fig = plotJustTheSpectra(scrsz, style, lambdaSimSPD, SPD, rea2005, criterionIndex, linLog)

    % Brainard and Thapan Data
    % [wavelength  response BrainarsdYES]
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
            
    %% FIGURE
    fig = figure('Color', 'w');
    
        if strcmp(linLog, 'lin') == 1
            set(fig, 'Position', [0.01*scrsz(3) 0.21*scrsz(4) 0.41*scrsz(3) 0.7*scrsz(4)])
        elseif strcmp(linLog, 'log') == 1
            set(fig, 'Position', [0.44*scrsz(3) 0.21*scrsz(4) 0.41*scrsz(3) 0.7*scrsz(4)])
        else
           error('linLog definition incorrect'); 
        end

        rows = 3;
        cols = 1;

        ocuMediaStrLabels = {'Aphakic', '25yr.', '65yr.'};

        %% RELATIVE
        jj = 1;
        indSp = jj;
        sp(indSp) = subplot(rows,cols,indSp);

            % Plot the melatonin suppression efficiency
            titleStr = 'Rea et al. (2005,2011)';
            yStr = sprintf('%s\n%s', 'Normalized', '1 / criterion response');
            p1 = plotTheMelSpectra(lambdaSimSPD, rea2005.ICLAeffNorm, style, titleStr, yStr, criterionIndex, linLog);      
            
                
                % Scaler, from Najjar et al.
                atWavelength = 460;
                wavelIndex = find(lambdaSimSPD == atWavelength);
                
                if ~isempty(wavelIndex)                   
                    scaler = rea2005.ICLAeffNorm{3}(wavelIndex);
                else
                    warning('Unscaled Najjar points, too coarse wavelength definition?')
                    scaler = 1;
                end

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

                % From Andrew's implementation
                load('ReaCurve.mat')
                if strcmp(linLog, 'lin') == 1
                    p1b(3) = plot(wave,ICLAeff,'b--');    
                elseif strcmp(linLog, 'log') == 1
                    p1b(3) = plot(wave,log10(ICLAeff),'b--'); 
                else
                   error('linLog definition incorrect');
                end                                            
                
                % Najjar
                if strcmp(linLog, 'lin') == 1
                    p1b(4) = errorbar(najjar(:,1), scaler*najjar(:,2), najjar(:,3), 'ko');
                elseif strcmp(linLog, 'log') == 1
                    logError_LO = log10(scaler*najjar(:,2) - najjar(:,3));
                    logError_HI = log10(scaler*najjar(:,2) + najjar(:,3));
                    p1b(4) = errorbar(najjar(:,1), log10(scaler*najjar(:,2)), logError_LO, logError_HI, 'ko');
                else
                    error('linLog definition incorrect');
                end
                set(p1b(4), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')

                hold off  
                leg(1) = legend('Aphakic', '25yr.', '65yr.', 'Brainard', 'Thapan', 'Orig. Curve', 'Najjar');
                legend('boxoff')
            

        %% REA ET AL. (2005,2011): ABSOLUTE
        jj = 2;
        indSp = jj;
        sp(indSp) = subplot(rows,cols,indSp);

            % Plot the melatonin suppression efficiency
            titleStr = 'Rea et al. (2005,2011)';
            yStr = sprintf('%s\n%s', 'Normalized to 25 yr Std. Observer', '1 / criterion response');

            % normalize the components before plotting
            norm1 = rea2005.normFactor{1} / rea2005.normFactor{2}
            rea2005.ICLAeffNorm{1}(:) = norm1 * rea2005.ICLAeffNorm{1}(:);                
            norm3 = rea2005.normFactor{3} / rea2005.normFactor{2}
            rea2005.ICLAeffNorm{3}(:) = norm3 * rea2005.ICLAeffNorm{3}(:); 
            
                % Scaler, from Najjar et al.
                atWavelength = 460;
                wavelIndex = find(lambdaSimSPD == atWavelength);
                
                if ~isempty(wavelIndex)                   
                    scaler = rea2005.ICLAeffNorm{3}(wavelIndex);
                else
                    warning('Unscaled Najjar points, too coarse wavelength definition?')
                    scaler = 1;
                end

            % plot
            p2 = plotTheMelSpectra(lambdaSimSPD, rea2005.ICLAeffNorm, style, titleStr, yStr, criterionIndex, linLog);     

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
                    
               
                if strcmp(linLog, 'lin') == 1
                    p2b(3) = errorbar(najjar(:,1), scaler*najjar(:,2), najjar(:,3), 'ko');
                elseif strcmp(linLog, 'log') == 1
                    logError_LO = log10(scaler*najjar(:,2) - najjar(:,3));
                    logError_HI = log10(scaler*najjar(:,2) + najjar(:,3));
                    p2b(3) = errorbar(najjar(:,1), log10(scaler*najjar(:,2)), logError_LO, logError_HI, 'ko');
                else
                    error('linLog definition incorrect');
                end
                set(p2b(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
                
                    leg(2) = legend('Aphakic', '25yr.', '65yr.', 'Brainard', 'Thapan', 'Najjar');
                    legend('boxoff')

        %% REA ET AL. (2005,2011): ABSOLUTE, long-wavelength normalization
        jj = 3;
        indSp = jj;
        sp(indSp) = subplot(rows,cols,indSp);

            % Plot the melatonin suppression efficiency
            titleStr = 'Rea et al. (2005,2011)';
            yStr = sprintf('%s\n%s', 'Normalized to 25 yr Std. Observer', '1 / criterion response');

                       
                % Scaler, from Najjar et al.
                atWavelength = 530;
                wavelIndex = find(lambdaSimSPD == atWavelength);
                
                if ~isempty(wavelIndex)                   
                    scaler = rea2005.ICLAeffNorm{3}(wavelIndex);
                else
                    warning('Unscaled Najjar points, too coarse wavelength definition?')
                    scaler = 1;
                end

            % plot
            p2 = plotTheMelSpectra(lambdaSimSPD, rea2005.ICLAeffNorm, style, titleStr, yStr, linLog);     

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
                    
                
                if strcmp(linLog, 'lin') == 1
                    p2b(3) = errorbar(najjar(:,1), scaler*najjar(:,2), najjar(:,3), 'ko');
                elseif strcmp(linLog, 'log') == 1
                    logError_LO = log10(scaler*najjar(:,2) - najjar(:,3));
                    logError_HI = log10(scaler*najjar(:,2) + najjar(:,3));
                    p2b(3) = errorbar(najjar(:,1), log10(scaler*najjar(:,2)), logError_LO, logError_HI, 'ko');
                else
                    error('linLog definition incorrect');
                end
                set(p2b(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
                
                    leg(2) = legend('Aphakic', '25yr.', '65yr.', 'Brainard', 'Thapan', 'Najjar');
                    legend('boxoff')

                    
                % oldCurve = norm3 * rea2005.ICLAeffNorm{3}(:);
                % dlmwrite('oldMelatoninCurve.txt',[lambda oldCurve])

            
    %% style

        % subplots (gca)
        set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)   
        set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)        
        set(sp, 'XLim', [400 620])
        if strcmp(linLog, 'lin') == 1
            set(sp, 'YLim', [0 1.2])
        elseif strcmp(linLog, 'log') == 1
            set(sp, 'YLim', [-2.5 0.6])
            set(leg, 'Location', 'South')
        end

        