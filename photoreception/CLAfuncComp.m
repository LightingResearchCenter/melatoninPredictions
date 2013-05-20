% Calculates CLA from SPD, v. postBerlinCorrMelanopsin_02Oct2012
function [CLA, Mel, SWS, Vl, Cones, Rod, CLAComp, oppFlag] = CLAfuncComp(spd, S_cornea, S_retina, spectralCrossover, modelType, a_cone, varargin)

    % Revised May 11, 2011
    % Calculates the circadian stimulus for the given spd
    % spd is assumed to be in units of W/m^2    
    % CS is scaled to be equal to lux for illuminanct A (2856 K)
    
    % spd is two column matrix with wavelength (nm) in the first column
    % and spectral irradiance (W/(m^2 nm) in the second column
    % OR spd is a column vector and start, end and increment wavelength values
    % are specified as additional arguements (e.g. f(spd,400,700,10))
    
    % Revised April, 2013
    % Modified by Petteri Teikari to work with lightLab library/functions
    % and to return also the spectral sensitivity curve ("CLASpec"), and
    % the basis vectors ("CLAComp")
    
    % For quick'n'dirty testing
    verifyModel = 1; % now the basis vectors are as in the original model
    
                     % when this is 1, then the spectral sensitivities are
                     % corrected for ocular media transmittance
                     
                     % when 0, the sensitivities are for retina, and then
                     % the ocular media is used to filter the light SPD
                     % making easier to simulate the effect of aged ocular
                     % media. Note in that case the V(l) and V(l)' are
                     % uncorrected for ocular media (25 yr standard
                     % observer)
                     
    %% Check input   
    
        % for quick'n'dirty development
        if nargin == 1
            load('specSensitivities.mat')
            wavelength_spd = (380:1:780)';
        else
            % save('specSensitivities.mat', 'S_cornea', 'S_retina')
        end    

        % Original checks from Andy
        if nargin > 1
            if isempty(varargin)
                [~, columns] = size(spd);
                if columns > 2
                    error('Not column oriented data. Try transposing spd');
                end
                wavelength_spd = (380:1:780)';
                spd = spd;
            else
                startw = varargin{1};
                endw = varargin{2};
                incrementw = varargin{3};
                wavelength_spd = (startw:incrementw:endw)';
                [rows columns] = size(spd);
                if columns > 1
                    error('Detected multiple columns of data. Try transposing spd');
                end
            end
        end
    
    %% Assign input variables in S to match the old notation    
        if verifyModel == 1        
            load('spectralSensitivitiesRea2011.mat')
            % here the original spectral sensitivities, for speeding this up
            
        else
            Vlambda = S_retina.vLambda;
                Vlambda(isnan(Vlambda)) = 0; % NaN -> 0
            Vprime = S_retina.vLambdaPrime;
            Scone = S_retina.SWS.data;
            M = S_retina.Melanop.data;
                % [~, Macula] = createMacularTemplate([380 780], 1); % max at 1 OD 
            % Macula = Macula / 2; % for peak density of 0.5 OD            
            %     macularTi = 10.^(-Macula); % transmittance
            %     save('macularTemplate.mat', 'macularTi')
            load('macularTemplate.mat')

            Scone = Scone ./ macularTi; % Correct for macular pigment
            Scone = Scone / (max(Scone)); % re-normalize
            Vlambda = Vlambda ./ macularTi; % Correct for macular pigment
            Vlambda = Vlambda / max(Vlambda); % re-normalize
        end        
        
    
    %% Define parameters
    
        rodSat = 35000; % Scotopic Trolands
        % retinalE = [1 3 10 30 100 300 1000 3000 10000 30000 100000];
        % pupilDiam = [7.1 7 6.9 6.8 6.7 6.5 6.3 5.65 5 3.65 2.3];
        % diam = interp1(retinalE, pupilDiam, rodSat, 'linear');
        diam = 3.5536; % mm at rodSat = 35000
        rodSat = rodSat / (diam^2 / 4*pi) * (pi/1700);

        % Define constants (Lighing Res. Tech 2011 notations)
        a1      = 1;        % 0.285
        b1      = 0.0;      % 0.01  
        a_by    = 0.6201;   % 0.2, a2
        b2      = 0.0;      % 0.001
        k       = 0.2616;   % 0.31
        a_rod   = 3.2347;   % 0.72, a3
        
        a_cone; % as input argument
        multip  = 1622.5;

    %% Define the basis vectors
            
        % Template vectors
        SCone_vec = Scone;
        VlambdaVec = Vlambda;
        MVec = M;
        VPrimeVec = Vprime;
        onesVec = ones(length(wavelength_spd),1);
        
        % Allow the intermediate results to be displayed when first
        % calculating the different component vectors rather than
        % integrating right away [trapz()]
        
            % The term for MELANOPSIN contribution
            CS1_Melanopsin = (a1 * (MVec.*spd));

            % The term for CONE contribution
            CS2_SCone = SCone_vec .*spd;
            CS2_LMCones = (k * (VlambdaVec.*spd));  
            CS2_Cones = a_by * ((CS2_SCone.*spd) - (CS2_LMCones.*spd));

            % The term defining ROD contribution
            Rod = a_rod * (1 - exp(-1 * ((VPrimeVec.*spd) / rodSat))); % for 'sharp'   
            rod_smoothAdd = onesVec - exp(40*(CS2_Cones));
            if strcmp(spectralCrossover, 'sharp')
                % add later this intermediate?                
                
            elseif strcmp(spectralCrossover, 'smooth')                
                Rod = Rod .* rod_smoothAdd;
            else
               error('error in defining spectralCrossover, maybe a typo on your code?') 
            end
    
        %% SPECTRAL OPPONENT RESPONSE        
        % trapz(wavelength_spd,CS2_SCone - CS2_LMCones)
        if trapz(wavelength_spd,CS2_SCone - CS2_LMCones) >= 0                      
                
            % Obtain scalars of different photoreceptor vectors 
            CS1scalar_Melanopsin = trapz(wavelength_spd,CS1_Melanopsin) - b1;
            
            CS2scalar_SCone = a_by * (trapz(wavelength_spd,Scone.*spd));
            CS2scalar_LMCones = a_by * (-k*trapz(wavelength_spd,Vlambda.*spd));
            
            CS2scalar_Cones = (CS2scalar_SCone + CS2scalar_LMCones) - b2; % trapz(wavelength_spd,CS2_Cones);
            
                % should be the same as CS2scalar_SCone + CS2scalar_LMCones
                % check = CS2scalar_Cones - (CS2scalar_SCone - CS2scalar_LMCones); % should be zero
                
            % The ROD scalar, with possibility to vary the spectral
            % crossover
            if strcmp(spectralCrossover, 'sharp') 
                RodScalar = -1*(a_rod*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)));
                
            elseif strcmp(spectralCrossover, 'smooth')
                RodScalar = -1*(a_rod*(1-exp(-trapz(wavelength_spd,Vprime.*spd)/rodSat)) ...
                    * exp(-40*(CS2scalar_Cones)));
            else
                error('error in defining spectralCrossover, maybe a typo on your code?') 
            end
            
                % Do the limit checks for the scalars
                if CS1scalar_Melanopsin < b1
                    CS1scalar_Melanopsin = 0; % remove negative values that are below threshold set by constant b1.
                end

                if CS2scalar_Cones < 0
                    CS2scalar_Cones = 0; % This is the important diode operator, the (b-y) term cannot be less than zero
                end            
            
            % Sum of different scalars, i.e. the CS
            CS = CS1scalar_Melanopsin + CS2scalar_Cones + RodScalar;            
                
                % Do the limit check for the final scalar
                if CS < 0
                    CS = 0; % Rod inhibition cannot make the CS less than zero
                end

            % disp('(B-Y) > 0')
            oppFlag = 1;
            
            % Assign the vectors to structur to be returned
            CSComp.CS1_Melanopsin = CS1_Melanopsin;
            CSComp.CS2_SCone = CS2_SCone;
            CSComp.CS2_LMCones = CS2_LMCones;
            CSComp.CS2_Cones = CS2_Cones;
            CSComp.Rod = Rod;
                
        %% MELANOPSIN-ONLY RESPONSE
        else

            CS1scalar_Melanopsin = trapz(wavelength_spd,CS1_Melanopsin) - b1;
            CS = CS1scalar_Melanopsin;  
            
            if CS < 0
                CS = 0; % Negative values mean stimulus is below threshold set by constant b1
            end
            
            if strcmp(modelType, 'original')  
                % nothing needs to be done
                
            elseif strcmp(modelType, 'origWithCones')                
                CS_coneDrive = a_cone * trapz(wavelength_spd,Vlambda.*spd);
                CS = CS_coneDrive + CS1scalar_Melanopsin;
            elseif strcmp(modelType, 'placeholderForModificationsForTheModel')
                % "doSomething" -block
            end

            %disp('(B-Y) < 0')
            oppFlag = 0;
            
            % Assign the vectors to structur to be returned
            CSComp.CS1_Melanopsin = CS1_Melanopsin;
            
            % these are NaN then
            CSComp.CS2_SCone = NaN;
            CSComp.CS2_LMCones = NaN;
            CSComp.CS2_Cones = NaN;
            CSComp.Rod = NaN;
            CS2scalar_SCone = NaN;
            CS2scalar_LMCones = NaN;
            CS2scalar_Cones = NaN;
            RodScalar = NaN; 
            
        end        
       
    %% FINALLY SCALE the output
    
        % Scalar
        CLA = multip * CS;

            % Scalar components
            Mel = multip * CS1scalar_Melanopsin;
            SWS = multip * CS2scalar_SCone;
            Vl  = multip * CS2scalar_LMCones;
            Cones = multip * CS2scalar_Cones;
            Rod = multip * RodScalar;    
        
        % Vectors
        CLAComp.CS1_Melanopsin =  multip * CSComp.CS1_Melanopsin;
        CLAComp.CS2_SCone = multip * CSComp.CS2_SCone;
        CLAComp.CS2_LMCones = multip * CSComp.CS2_LMCones;
        CLAComp.CS2_Cones = multip * CSComp.CS2_Cones;
        CLAComp.Rod = multip * CSComp.Rod;
            
            
function [mp_density, mp_density_norm] = createMacularTemplate(lambdaLimits, lambda_res)

    % The model derived by by Walraven [1] as referred by Zagers and van
    % Norren (2004) [2]; and van de Kraats and van Norren (2008) [3].
    
    % INPUTS:
    %
    %   lambdaLimits  2-element vector where the first "lambda(1)" is the
    %                 minimum wavelength in the created template and second
    %                 value "lambda(2)" is the maximum wavelength of the
    %                 created vector. Wavelengths are given in nanometers
    %
    %   lambda_res    the spectral resolution of the created template in
    %                 nanometers
    %
    % OUTPUTS:
    %

    %
    % EXAMPLE OF USE:
    %
    %   macularPigmentOD = createMacularTemplate([380 780], 1)
    %
    % REFERENCES:
    %
    %   [1] P. L. Walraven, CIE Rep. TC 1-36, Draft 7
    %
    %   [2] Niels P. A. Zagers and Dirk van Norren, “Absorption of the eye 
    %       lens and macular pigment derived from the reflectance of cone photoreceptors,” 
    %       Journal of the Optical Society of America A 21, no. 12 (December 1, 2004): 2257-2268.
    %       http://dx.doi.org/10.1364/JOSAA.21.002257
    % 
    %   [3] Jan van de Kraats and Dirk van Norren, “Directional and nondirectional 
    %       spectral reflection from the human fovea,” 
    %       Journal of Biomedical Optics 13, no. 2 (March 0, 2008): 024010-13.
    %       http://dx.doi.org/10.1117/1.2899151
    
    % create the wavelength vector
    lambda  = (linspace(lambdaLimits(1), lambdaLimits(2), (((lambdaLimits(2)-lambdaLimits(1)) / lambda_res) + 1)) )';
    
    % create the scalar vector needed for vectorized version of the
    % template creation
    ones436 = 436 * ones(length(lambda(:,1)), 1);
    ones480 = 480 * ones(length(lambda(:,1)), 1); 
    ones458 = 458 * ones(length(lambda(:,1)), 1);
    ones457 = 457 * ones(length(lambda(:,1)), 1);    
    
    % create the template using the given equation [1] and for the
    % wavelengths specified just above  
    opticalDensity =  (0.32 .* exp(-0.0012 .* ((ones436 - lambda) .^ 2))) ... 
                   + (0.32 .* exp(-0.0012 .* ((ones480 - lambda) .^ 2)))...
                   - (0.123 .* exp(-0.0012 .* ((ones458 - lambda) .^ 2))) ...
                   + (0.12042 .* exp(-0.006 .* ((ones457 - lambda) .^ 2)));

   % assign the results to output variable   
   mp_density = opticalDensity;
   
   % normalize to unity at 460 nm                
   mp_density_norm = mp_density / max(mp_density);   
   
        % if you want to normalize to unity at 460 nm just multiply the created
        % template with a scalar f_norm which was set to 1/0.35 by Zagers and
        % van Norren [2].   
            % f_norm = 1 / 0.35;
            % template(:,1) = template(:,1);
            % templatem(:,2) = f_norm * template(:,2);