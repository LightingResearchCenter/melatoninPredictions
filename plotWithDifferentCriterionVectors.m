function fig = plotWithDifferentCriterionVectors(scrsz, style, lambda, SPD, rea2005, criteriaVector, linLog)

    fig = figure('Color', 'w');
        set(fig, 'Position', [0.01*scrsz(3) 0.51*scrsz(4) 0.98*scrsz(3) 0.4*scrsz(4)])
        
        % subplot layout
        rows = 1;
        cols = 3;
        
    % normalize the components before plotting
    norm1 = rea2005.normFactor{1} / rea2005.normFactor{2};
    rea2005.ICLAeff{1}(:) = norm1 * rea2005.ICLAeff{1}(:);                
    norm3 = rea2005.normFactor{3} / rea2005.normFactor{2};
    rea2005.ICLAeff{3}(:) = norm3 * rea2005.ICLAeff{3}(:);

    yStr = {'Aphakic'; '25yr'; '65yr'};
    
    % Absolute values
    for ij = 1 : length(rea2005.ICLAeff)
    
        sp(ij) = subplot(rows,cols,ij);            
        
            legendStr = cell(length(criteriaVector),1);
            
            for ik = 1 : length(criteriaVector)
            
                if strcmp(linLog, 'lin') == 1
                    hold on
                    p(ij,ik) = plot(lambda, rea2005.ICLAeff{ij}(:,ik), 'b');
                    legendStr{ik} = ['C=', num2str(criteriaVector(ik))];
                                     
                elseif strcmp(linLog, 'log') == 1
                    hold on
                    p(ij,ik) = plot(lambda, log10(rea2005.ICLAeff{ij}(:,ik)), 'b');
                    legendStr{ik} = ['C=', num2str(criteriaVector(ik))];
                else
                    error('linLog definition incorrect');
                end   
                
            end    
            
            lab(ij,1) = xlabel('Wavelength [nm]');     
                lab(ij,2) = ylabel(yStr{ij});
            
            leg(ij) = legend(legendStr);
                    legend('boxoff')    
                    
            if ij == 2
                titleStr = sprintf('%s\n%s', 'Melatonin suppression', 'as a function of criterion response threshold');
                tit = title(titleStr);  
            end
            
    end
    
    % style
    set(p(:,2), 'Color', 'k')
    set(p(:,3), 'Color', 'r')
    set(p(:,4), 'Color', [0.3 0.8 0])
    set(p(:,5), 'Color', [0 0.5 0.75])
    
    set(sp, 'XLim', [400 650])
    
    set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)   
    set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)  
    set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')  
    set(tit, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')

    
    
    