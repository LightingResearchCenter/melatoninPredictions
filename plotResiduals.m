function plotResiduals(brainardStats, thapanStats, najjarStats, titleStr, handles)

    if nargin == 0
        load tmmmmmp.mat
    else
        save tmmmmmp.mat
    end

    residMatrixBrainard = zeros(length(brainardStats{1}.residuals.x), length(brainardStats)+1);
    residMatrixThapan = zeros(length(thapanStats{1}.residuals.x), length(brainardStats)+1);
    residMatrixNajjar = zeros(length(najjarStats{1}.residuals.x), length(brainardStats)+1);

    for i = 1 : length(brainardStats)
              
        % brainard        
        x = brainardStats{i}.residuals.x;
        y = brainardStats{i}.residuals.y;
        
        residMatrixBrainard(:,1) = x;
        residMatrixBrainard(:,i+1) = y;       
        
        % thapan        
        x = thapanStats{i}.residuals.x;
        y = thapanStats{i}.residuals.y;
        residMatrixThapan(:,1) = x;
        residMatrixThapan(:,i+1) = y;        
        
        % najjar
        x = najjarStats{i}.residuals.x;
        y = najjarStats{i}.residuals.y;
        residMatrixNajjar(:,1) = x;
        residMatrixNajjar(:,i+1) = y;                
        
    end
    
    numberOfPoints = length(brainardStats{i}.residuals.x) + length(thapanStats{i}.residuals.x) + length(najjarStats{i}.residuals.x);    
    residMatrixAll = zeros(numberOfPoints, length(brainardStats)+1);
    
    for i = 1 : length(brainardStats)        
        
        % brainard        
        x = brainardStats{i}.residuals.x;
        y = brainardStats{i}.residuals.y;
        ind1(1) = 1;
        ind2(1) = length(x);
        residMatrixAll(ind1(1):ind2(1),1) = x;
        residMatrixAll(ind1(1):ind2(1),i+1) = y;
                
        % thapan        
        x = thapanStats{i}.residuals.x;
        y = thapanStats{i}.residuals.y;
        ind1(2) = ind2(1) + 1;
        ind2(2) = ind1(2) + length(x) - 1;
        residMatrixAll(ind1(2):ind2(2),1) = x;
        residMatrixAll(ind1(2):ind2(2),i+1) = y;
       
        % najjar
        x = najjarStats{i}.residuals.x;
        y = najjarStats{i}.residuals.y;
        ind1(3) = ind2(2) + 1;
        ind2(3) = ind1(3) + length(x) - 1;
        residMatrixAll(ind1(3):ind2(3),1) = x
        residMatrixAll(ind1(3):ind2(3),i+1) = y;               
        
    end
      
    % sort
    [residMatrixAllSorted,index] = sortrows(residMatrixAll(ind1(1):ind2(2),:))
    
    sp(1) = subplot(2,1,1)
        pH = plot(residMatrixAllSorted(:,1), abs(residMatrixAllSorted(:,2:end)), '-o');
            legend(titleStr)
            legend('boxoff')

            set(pH, 'MarkerSize', 4)
            xlim([400 650])
            title('25yr : Brainard/Thapan')
        
    sp(2) = subplot(2,1,2)
        pH = plot(residMatrixAll(ind1(3):ind2(3),1), abs(residMatrixAll(ind1(3):ind2(3),2:end)), '-o');
            legend(titleStr)
            legend('boxoff')

            set(pH, 'MarkerSize', 4)
            xlim([400 650])
            title('65yr : Najjar')
        
        
    % FIRST LOOP
    %{
    hold on
    pH = plot(residMatrixBrainard(:,1), abs(residMatrixBrainard(:,2:end)), '-bo', ...
              residMatrixThapan(:,1), abs(residMatrixThapan(:,2:end)), '-ro', ...
              residMatrixNajjar(:,1), abs(residMatrixNajjar(:,2:end)), '-go')
    hold off    
        
        lineColors = [0 0 0; 0 1 0; 1 1 0; 1 0 1];
    
        % set line colors
        for i = 1 : length(brainardStats)
            % set(pH([1 5 9 13]), 'Color', lineColors(2,:))
        end
        
        % set marker colors
        set(pH(1:4), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
        set(pH(5:8), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
        set(pH(9:12), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
        
        legend(titleStr)
    %}
        
        