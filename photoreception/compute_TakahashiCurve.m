function [pupil, melSuppr] = compute_TakahashiCurve(lambdaSimSPD, SPD, irrad, templates)

    nrOfIterations = length(SPD) * length(lambdaSimSPD);
    for j = 1 : length(SPD) % no of ocular media models

        for i = 1 : length(lambdaSimSPD) % no of light SPD
            
            k = 20;           
            loopIndex = ((i-j) * length(lambdaSimSPD)) + i;
            disp(['', num2str(100*(loopIndex/nrOfIterations),'%3.1f'), '%'])
            
            [pupilRaw, melSupprRaw] = create_TakahashiCurve(SPD{j}{i}{k}, templates);
            
            pupil{j}(i) = pupilRaw.Area;
            melSuppr{j}(i) = melSupprRaw.dilated;
            
        end
        
    end