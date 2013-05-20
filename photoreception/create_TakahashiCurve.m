function [pupil, melSuppr] = create_TakahashiCurve(SPD, templatePath)

    %     Takahashi, Yoshika, Tetsuo Katsuura, Yoshihiro Shimomura, and Koichi Iwanaga. 
    %     “Prediction Model of Light-induced Melatonin Suppression.” 
    %     Journal of Light & Visual Environment 35, no. 2 (2011): 123–135. 
    %     http://dx.doi.org/10.2150/jlve.35.123

    % Implementation by 
    % Petteri Teikari, 2013, LRC, Troy, NY, USA
    % teikap@rpi.edu

    if nargin == 0
        templatePath = '/home/petteri/Dropbox/Matlab Code/In Development/lightLab/database/Templates';
        dataRaw = importdata(fullfile(templatePath, 'takahashi2011_CieA.txt'));
            wavelength = dataRaw.data(:,1);
            SPD = dataRaw.data(:,2); % Irradiance [10-6 W cm-2]
    end
    
    %% import the data   

        dataRaw = importdata(fullfile(templatePath, 'takahashi2011_calcSheet.txt'), '\t', 1);
            wavelength = dataRaw.data(:,1);
            photonE = dataRaw.data(:,2); % Energy of one photon [J]
            
            % this need to be overwritten if you use some other light
            % source
            photonDensity = dataRaw.data(:,3); % Photon density [1012 photons cm-2 s-1]
            illumPhotopic = dataRaw.data(:,4); % Irradiance [W/m2] * V
            illumPhotopic = dataRaw.data(:,5); % Irradiance [W/m2] * V'
            illumMelanopic = dataRaw.data(:,6); % Irradiance [W/m2] * M(l)
            illumIPRGC = dataRaw.data(:,7); % Irradiance [W/m2] * ipRGC(l)

        dataRaw = importdata(fullfile(templatePath, 'takahashi2011_basicData.txt'), '\t', 1);
            
            wavelength = dataRaw.data(:,1);
            vlambda = dataRaw.data(:,2); 
            vLambdaScotopic = dataRaw.data(:,3); 
            mLambda = dataRaw.data(:,4); 
            ipLambda = dataRaw.data(:,5); 
    
    %% Calculate dot products, weigh SPD with sensitivity functions
        
        [pupil, melSuppr] = takahashi_weighSPD(SPD, photonE, vlambda, vLambdaScotopic, mLambda, ipLambda);    
    
        
    function [pupil, melSuppr] = takahashi_weighSPD(SPD, photonE, vlambda, vLambdaScotopic, mLambda, ipLambda)
    
        % Photon density [1012 photons cm-2 s-1]
        photonDensity.vector = (SPD * 10^-6) ./ photonE / 10^12;
        photonDensity.scalar = sum(photonDensity.vector);
        
        % Irradiance [W/m2] * V  
        illumPhotopic.vector = (SPD * 10^-6 / 10^-4) .* vlambda;
        illumPhotopic.scalar = 683 * sum(illumPhotopic.vector);
        
        % Irradiance [W/m2] * V'
        illumScotopic.vector = (SPD * 10^-6 / 10^-4) .* vLambdaScotopic;
        illumScotopic.scalar = 1700 * sum(illumPhotopic.vector);
        
        % Irradiance [W/m2] * M(l)
        illumMelanopic.vector = (SPD * 10^-6 / 10^-4) .* mLambda;
        illumMelanopic.scalar = 5450 * sum(illumMelanopic.vector);
        
         % Irradiance [W/m2] * ipRGC(l)
        illumIPRGC.vector = (SPD * 10^-6 / 10^-4) .* ipLambda;
        illumIPRGC.scalar = 3330 * sum(illumIPRGC.vector);
        
        ipRGC = illumIPRGC.scalar;
        
    %% Calculate pupil size
    
        pupil.Diam = 7.75 - (5.75 * ( ((4135*illumScotopic.scalar /1692 ) ^0.41) / ( (4135*illumScotopic.scalar /1692 )^0.41+2))); % not correct
        pupil.Area = (pupil.Diam/2)^2 * pi();
        
    %% MELATONIN SUPPRESSION
    
        alpha = ((0-1) / ((1 + (illumIPRGC.scalar / 845) ^ 2.26)) + 1);
    
        maxSuppr = 66.9;
        slope = 1.27;
        x50 = 131;
        
        melSuppr.dilated = (0 - maxSuppr) / (1 + ((illumMelanopic.scalar/x50)^slope)) + maxSuppr;
        melSuppr.wPupil = (0-maxSuppr) / (1 + ( (illumMelanopic.scalar*pupil.Area) / (x50 * (7.19/2)^2 * pi())) ^slope ) + maxSuppr;
        melSuppr.aNormal = alpha * melSuppr.wPupil;
        
  
      