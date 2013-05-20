function rea2005 = computeReaParameters(lambdaSimSPD, SPD, irrad, criteriaVector, S_cornea, S_retina, spectralCrossover, modelType, a_cone, path)

    % you could probably vectorize this for better computational efficiency

    % for fine spectral and/or irradiance resolution, computation
    % might take some time so it is nice maybe to keep track of the
    % progress
    nrOfIterations = 3 * length(lambdaSimSPD) * length(irrad); 
    cd(path.library)    
    
    % preallocate
    CLA = cell(length(SPD),1);
    Mel = cell(length(SPD),1);
    SWS = cell(length(SPD),1);
    Vl = cell(length(SPD),1);
    Cones = cell(length(SPD),1);
    Rod = cell(length(SPD),1);
    CLAComp = cell(length(SPD),1);
    oppFlag = zeros(length(SPD),length(irrad));
    
    for j = 1 : length(SPD) % no of ocular media models

        for i = 1 : length(lambdaSimSPD) % no of light SPD            
            
            for k = 1 : length(irrad) % criterion/irradiance
                
                % Compute Rea et al. 2005 response
                cd(path.photoreception)                             
                    [CLA{j}(k), Mel{j}(k), SWS{j}(k), Vl{j}(k), Cones{j}(k), Rod{j}(k), CLAComp{j}{k}, oppFlag(j,k)] ...
                        = CLAfuncComp(SPD{j}{i}{k}, S_cornea, S_retina, spectralCrossover, modelType, a_cone);
               
            end
            
            % "progress bar"
            loopIndex = ((j-1) * length(lambdaSimSPD) * length(irrad)) + ((i-1) * length(irrad)) + k;
            disp(['', num2str(100*(loopIndex/nrOfIterations),'%3.1f'), '%'])

            % find the irradiance that gives the criterion response     
            for kj = 1 : length(criteriaVector)
                
                try             

                    rea2005.ICLA{j}(i,kj) = interp1(CLA{j},irrad, criteriaVector(kj),'linear');          

                    % Feed now this irradiance value back to the
                    % CLAfuncComp() to obtain the corresponding
                    % component values to this irradiance

                        % Create the modified light SPD
                        SPD_atCriterion = SPD{j}{i}{k} / sum(SPD{j}{i}{k}); % normalize to one
                        SPD_atCriterion = SPD_atCriterion * rea2005.ICLA{j}(i,kj); % correct the shape to correct irradiance                              

                            % a = rea2005.ICLA{j}(i);
                            % plot(SPD_atCriterion); title(num2str(a))
                            % drawnow; pause(0.01)

                        % Call the function
                        [CLA2{j}(i), Mel2{j}(i), SWS2{j}(i), Vl2{j}(i), Cones2{j}(i), Rod2{j}(i), CLAComp2{j}, oppFlag2(j)] ...
                            = CLAfuncComp(SPD_atCriterion, S_cornea, S_retina, spectralCrossover, modelType, a_cone);

                        % Assign to a structure                              
                        rea2005.Mel{j}(i,kj) = Mel2{j}(i);
                        rea2005.SWS{j}(i,kj) = SWS2{j}(i);
                        rea2005.Vl{j}(i,kj) = Vl2{j}(i);
                        rea2005.Cones{j}(i,kj) = Cones2{j}(i);
                        rea2005.Rod{j}(i,kj) = Rod2{j}(i);                                                      

                        % Sum together
                        rea2005.compSum{j}(i,kj) = nansum([rea2005.Mel{j}(i,kj) rea2005.Cones{j}(i,kj) rea2005.Rod{j}(i,kj)]);       

                            if rea2005.compSum{j}(i,kj) == 0
                               rea2005.compSum{j}(i,kj) = NaN; 
                            end             

                catch err % the interpolation will fail if all the CLA values are zero           

                    % err.identifier

                    if (strcmp(err.identifier,'MATLAB:griddedInterpolant:NonMonotonic1DErrId'))
                        warning(err.identifier,['all CLA values are zeroes, thus the ICLA(i)=NaN, loop indices = ', num2str([j i k loopIndex nrOfIterations])])
                        rea2005.ICLA{j}(i,kj) = NaN;
                        rea2005.Mel{j}(i,kj)  = NaN;
                        rea2005.SWS{j}(i,kj)  = NaN;
                        rea2005.Vl{j}(i,kj)   = NaN;
                        rea2005.Cones{j}(i,kj) = NaN;
                        rea2005.Rod{j}(i,kj)  = NaN;
                        rea2005.compSum{j}(i,kj) = NaN;

                    end 
                end                       
            end
        end
        
        % go through the criteria vector
        for kj = 1 : length(criteriaVector)
        
            % CLA

                % efficiency  = 1 / (criterion response)                
                rea2005.ICLAeff{j}(:,kj) = 1 ./ rea2005.ICLA{j}(:,kj);                                        

                % normalize to a maximum of 1
                rea2005.normFactor{j}(kj) = max(rea2005.ICLAeff{j}(:,kj));                    
                rea2005.ICLAeffNorm{j}(:,kj) = rea2005.ICLAeff{j}(:,kj) / rea2005.normFactor{j}(kj);

            % CLA Components

                % efficiency  = 1 / (criterion response)                           
                rea2005.IMel{j}(:,kj) = 1 ./ rea2005.Mel{j}(:,kj);
                rea2005.ISWS{j}(:,kj) = 1 ./ rea2005.SWS{j}(:,kj);
                rea2005.IVl{j}(:,kj) = 1 ./ rea2005.Vl{j}(:,kj);
                rea2005.ICones{j}(:,kj) = 1 ./ rea2005.Cones{j}(:,kj);
                rea2005.IRod{j}(:,kj) = 1 ./ rea2005.Rod{j}(:,kj);

                % normalize when plotting, below plotTheComponents()
                
        end

    end
    cd(path.mainCode)  