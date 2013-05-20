function p = plotTheMelSpectra(lambda, supprEfficiency, style, titleStr, yStr, criterionIndex, linLog)
 

    % plot
    hold on
    if strcmp(linLog, 'lin') == 1
        p(1:3) = plot(lambda, supprEfficiency{1}(:,criterionIndex), 'b', ...
                      lambda, supprEfficiency{2}(:,criterionIndex), 'k', ...
                      lambda, supprEfficiency{3}(:,criterionIndex), 'y');  
    elseif strcmp(linLog, 'log') == 1
        p(1:3) = plot(lambda, log10(supprEfficiency{1}(:,criterionIndex)), 'b', ...
                      lambda, log10(supprEfficiency{2}(:,criterionIndex)), 'k', ...
                      lambda, log10(supprEfficiency{3}(:,criterionIndex)), 'y');  
    else
        error('linLog definition incorrect');
    end
    
    % save to disk, tabular data
    % NOTE! filename is now fixed
    dlmwrite('reaEtAl_380to780nm_1nm_hbw10nm.txt', [lambda supprEfficiency{2}(:,criterionIndex)], 'delimiter', '\t')
              
                

    lab(1) = xlabel('Wavelength [nm]');     
    lab(2) = ylabel(yStr);
    tit = title(titleStr);    
    
    leg = legend('Aphakic', '25yr.', '65yr.'); 
    legend('boxoff');
    
    % line plots            
    set(p(:,1), 'Color', style.colorPlot(1,:))
    set(p(:,2), 'Color', style.colorPlot(2,:))
    set(p(:,3), 'Color', style.colorPlot(3,:))

    % legend % text
    set(leg, 'Location', 'NorthEast', 'FontName', style.fontName, 'FontSize', style.fontBaseSize-1)        
    set(lab, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')    
    set(tit, 'FontName', style.fontName, 'FontSize', style.fontBaseSize+1, 'FontWeight', 'bold')
    