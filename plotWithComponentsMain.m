function fig = plotWithComponentsMain(scrsz, style, lambdaSimSPD, SPD, rea2005, criterionIndex, linLog)

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

    %% FIGURE
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.01*scrsz(3) 0.21*scrsz(4) 0.98*scrsz(3) 0.7*scrsz(4)])

        rows = 2;
        cols = 5;

        ocuMediaStrLabels = {'Aphakic', '25yr.', '65yr.'};

        %% RELATIVE
        jj = 1;
        indSp = ((jj-1)*cols)+1;
        sp(indSp) = subplot(rows,cols,[indSp indSp+1]);

            % Plot the melatonin suppression efficiency
            titleStr = 'Rea et al. (2005,2011)';
            yStr = sprintf('%s\n%s', 'Normalized', '1 / criterion response');
            p1 = plotTheMelSpectra(lambdaSimSPD, rea2005.ICLAeffNorm, style, titleStr, yStr, criterionIndex, linLog);      

                % Brainard & Thapan data                
                if strcmp(linLog, 'lin') == 1
                    p1b(1) = plot(waveBrainard, Brainard, 'bo');
                    p1b(2) = plot(waveThapan, Thapan, 'ro');
                elseif strcmp(linLog, 'log') == 1
                    p1b(1) = semilogy(waveBrainard, Brainard, 'bo');
                    p1b(2) = semilogy(waveThapan, Thapan, 'ro');
                else
                   error('linLog definition incorrect'); 
                end

                % From Andy's implementation
                load('ReaCurve.mat')
                if strcmp(linLog, 'lin') == 1
                    p1b(3) = plot(wave,ICLAeff,'b--');    
                elseif strcmp(linLog, 'log') == 1
                    p1b(3) = semilogy(wave,ICLAeff,'b--'); 
                else
                   error('linLog definition incorrect');
                end
                hold off              

                leg = legend('Aphakic', '25yr.', '65yr.', 'Brainard', 'Thapan', 'Orig. Curve');

            % Plot the components                
            for ij = 1 : length(SPD)      
                indSp = ((jj-1)*cols)+2+ij;                
                sp(indSp-1) = subplot(rows,cols,indSp);
                [pComp1(ij,:), legComp1(ij)] = plotTheComponents(lambdaSimSPD,rea2005,criterionIndex,ij,jj,0);
                titComp(jj,ij) = title(ocuMediaStrLabels{ij});
                if ij == 1
                    str = sprintf('%s\n%s', 'Normalized', 'criterion response');
                    lab2(1) = ylabel(str);
                end
            end

        %% REA ET AL. (2005,2011): ABSOLUTE
        jj = 2;
        indSp = ((jj-1)*cols)+1;
        sp(indSp-1) = subplot(rows,cols,[indSp indSp+1]);

            % Plot the melatonin suppression efficiency
            titleStr = 'Rea et al. (2005,2011)';
            yStr = sprintf('%s\n%s', 'Normalized to 25 yr Std. Observer', '1 / criterion response');

            % normalize the components before plotting
            norm1 = rea2005.normFactor{1} / rea2005.normFactor{2};
            rea2005.ICLAeffNorm{1}(:) = norm1 * rea2005.ICLAeffNorm{1}(:);                
            norm3 = rea2005.normFactor{3} / rea2005.normFactor{2};
            rea2005.ICLAeffNorm{3}(:) = norm3 * rea2005.ICLAeffNorm{3}(:); 

            % plot
            p2 = plotTheMelSpectra(lambdaSimSPD, rea2005.ICLAeffNorm, style, titleStr, yStr, criterionIndex, linLog);            

                % oldCurve = norm3 * rea2005.ICLAeffNorm{3}(:);
                % dlmwrite('oldMelatoninCurve.txt',[lambda oldCurve])

            % Plot the components
            for ij = 1 : length(SPD)      
                indSp = ((jj-1)*cols)+2+ij;
                sp(indSp-2) = subplot(rows,cols,indSp);
                [pComp2(ij,:), legComp2(ij)] = plotTheComponents(lambdaSimSPD,rea2005,criterionIndex,ij,jj,1);
                titComp(jj,ij) = title(ocuMediaStrLabels{ij});
                if ij == 1
                    str = sprintf('%s\n%s', 'Normalized to 25 yr Std. Observer', 'criterion response');
                    lab2(2) = ylabel(str);
                end
            end

    %% style

        % subplots (gca)
        set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)        
        set(sp, 'XLim', [400 600])

        set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)
        set([legComp1 legComp2], 'Location', 'SouthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)
        set([legComp1(3) legComp2(3)], 'Location', 'NorthEast')
        set(lab2, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')
        set(titComp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')