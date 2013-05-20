 function  [p, leg] = plotTheComponents(lambdaSimSPD, rea2005, criterionIndex, ij, jj, differFlag)            
        
        yMatrix = [rea2005.Mel{ij}(:,criterionIndex) rea2005.SWS{ij}(:,criterionIndex) rea2005.Vl{ij}(:,criterionIndex) ...
                   rea2005.Cones{ij}(:,criterionIndex) rea2005.Rod{ij}(:,criterionIndex) rea2005.compSum{ij}(:,criterionIndex)];    
           

        % yMatrix = [rea2005.IMel{ij}' rea2005.ISWS{ij}' rea2005.IVl{ij}' rea2005.ICones{ij}' rea2005.IRod{ij}'];
        % yMatrix(isinf(yMatrix)) = NaN;

    if differFlag == 1

        ind = 2; % 25yr
        normMatrix = [rea2005.Mel{ind}(:,criterionIndex) rea2005.SWS{ind}(:,criterionIndex) rea2005.Vl{ind}(:,criterionIndex) ...
                      rea2005.Cones{ind}(:,criterionIndex) rea2005.Rod{ind}(:,criterionIndex) rea2005.compSum{ij}(:,criterionIndex)];
        % normMatrix = [rea2005.IMel{ind}' rea2005.ISWS{ind}' rea2005.IVl{ind}' rea2005.ICones{ind}' rea2005.IRod{ind}'];
        % normMatrix(isinf(normMatrix)) = NaN;

        norm = nanmax(nanmax(normMatrix));

        % i2 = ij;
        % if ij == 3                
        %     i1 = 1; % aphakic
        % else
        %     i1 = ij + 1;
        % end            
        % yCol1 = [rea2005.IMel{i1}' rea2005.ISWS{i1}' rea2005.IVl{i1}' rea2005.ICones{i1}' rea2005.IRod{i1}'];
        % yCol2 = [rea2005.IMel{i2}' rea2005.ISWS{i2}' rea2005.IVl{i2}' rea2005.ICones{i2}' rea2005.IRod{i2}'];         
        % yMatrix = yCol1 ./ yCol2;            
        % pause         

    else     
        % a = 2;
        norm = nanmax(nanmax(yMatrix));
        % norm = 1;
    end
    yMatrix = yMatrix / norm;

    p = plot(lambdaSimSPD, yMatrix);
        legStr = {'Melanopsin'; 'S-cone'; 'V(\lambda)'; 'Cones'; 'Rod'; 'CLA'};
        leg = legend(legStr);
            legend('boxoff')
            box off


