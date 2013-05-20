function [raw, orig, irrad] = createMonochromaticLightSPDs(lambda, lambdaSimSPD, FWHM, photDensityCnst, irradMin, irradMax, irradRes, path)    

    xRes = lambda(2) - lambda(1);
    xLimits = [min(lambda) max(lambda)];   

    % test the model with different intensities  (W/m^2)
    % -4:0.2:3; % steps of log10(irradiance)
    irrad = (irradMin : irradRes : irradMax)';
    irrad = 10 .^ irrad;
    
    % pause
    
    cd(path.nomogram)
    for i = 1 : length(lambdaSimSPD)
        
        if length(FWHM) == 1 % FWHM the same for all wavelengths
            FWHM_loop = FWHM;
        else % and if the light have different hbws
            FWHM_loop = FWHM(i);
        end

        % generate the shape
        raw{i} = monochromaticLightAsGaussian(lambdaSimSPD(i), FWHM_loop, xRes, xLimits);

        % total irradiance
        for k = 1 : length(irrad)
            raw{i} = raw{i} / sum(raw{i}); % normalize to a radiant power of 1
            orig{i}{k} = irrad(k) * raw{i}; % weigh with the irradiance vector             
            
            % if you wanna do constant photon density        
            if photDensityCnst == 1
                cd(path.common)
                Q = convert_fromEnergyToQuanta(raw{i}, lambda);
                norm = max(Q);
                Q = Q / norm; % norm to one
                orig{i}{k} = convert_fromQuantaToEnergy(Q, lambda);
            end
            
        end 
        

    end
    cd(path.mainCode)