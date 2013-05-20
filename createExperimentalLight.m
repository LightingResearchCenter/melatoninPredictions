function SPD = createExperimentalLight(lambdaSimSPD, lambda, SPDraw, wavelengthInd, photonDensity)

    % input as photon density, output as irradiance, uW/cm^2
    %     
    %     SPDraw
    %     lambdaSimSPD
    %     wavelengthInd
    %     photonDensity    
    %     whos
    
    % Now scale the monochromatic spectra so that they have the desired
    % photon density
    SPD1 = SPDraw{1}{wavelengthInd(1)};
        multiplier1 = photonDensity(1) / sum(SPD1); % could use trapz as well
        SPD1 = multiplier1 * SPD1;
        sumNow1 = sum(SPD1); % debug
        
        if length(wavelengthInd) > 1
            SPD2 = SPDraw{1}{wavelengthInd(2)};
            multiplier2 = photonDensity(2) / sum(SPD2); % could use trapz as well
            SPD2 = multiplier2 * SPD2;
            sumNow2 = sum(SPD2); % debug
        else
            SPD2 = zeros(length(SPD1),1);
        end
    
    % Sum now these together (bichromatic stimulus)
    SPD_photonDensity = SPD1 + SPD2;
    photonDensitySum = sum(SPD_photonDensity);
    
    % and finally convert to irradiance (W/cm^2)
    SPD = convert_fromQuantaToEnergy(SPD_photonDensity, lambda);
    
    % Rea model requires W/m^2, thus
    SPD = SPD * 10^4;
    
    irradiance = sum(SPD);
    
    
    
    function E = convert_fromQuantaToEnergy(Q, lambda)

        
     % Inputs
     %      Q       - Photon density [ph/cm^2/sec]
     %      lambda  - Wavelength [nm]
     % Output
     %      E       - Irradiance [W/cm^2/sec (/nm)]
     
     h = 6.62606896 * 10^-34; % Planck's constant [J*s]
     c = 299792458; % Speed of light [m/s]  
     photonEnergy_vector = (h * c) ./ (lambda * 10^-9); % [J]
     %whos
     E = Q .* photonEnergy_vector;