function S = nomog_Govardovskii2000(peak_nm, bands, linLog, quantaE, density, xRes, xLims, optFilter)

    %% INPUTSs:
    % =======
    
        % peak_nm  - is the wavelength of peak sensitivity in nm. 
        % linLog   - output, LINEAR or LOG
        % quantaE  - output, in QUANTA or as ENERGY
        % density  - peak axial density for self-screening
        % xLims    - min/max wavelength for the created vector                
        % xRes     - the spectral resolution of the created nomogram [nm]      
        
            % Default density values from the literature
            % 0.50 in Tsujimura et al. (2010) for melanopsin,   http://dx.doi.org/10.1098/rspb.2010.0330
            % 0.27 in Lamb (1995) CONES, (diff. for S-cones!)   http://dx.doi.org/10.1016/0042-6989(95)00114-F 
            % 0.38 for M/L-cones, Stockman et al. (1999),       http://dx.doi.org/10.1016/S0042-6989(98)00225-9
            % 0.30 for S-cones, Stockman et al. (1999),         http://dx.doi.org/10.1016/S0042-6989(98)00225-9
            % 0.4 in Lamb (1995) RODS,                          http://dx.doi.org/10.1016/0042-6989(95)00114-F        

            
    %% OUTPUT:
    % =======
        % S - vector, continuous action spectrum
        
        
    %% CHECK for NUMBER of inputs    
    
        % default values
        def_bands = 'both';
        def_linLog  = 'lin';        
        def_quantaE = 'Q';
        def_density = 0.01;
        def_xRes    = 1;
        def_xLims   = [380 780];
        def_optFilter = [];

        if nargin == 0
            errordlg('Supply at least the peak wavelength')                
        elseif nargin == 1
            bands   = def_bands;
            linLog  = def_linLog;
            quantaE = def_quantaE;
            density = def_density;
            xRes    = def_xRes;
            xLims   = def_xLims;
            optFilter = def_optFilter;
        elseif nargin == 2        
            linLog  = def_linLog;
            quantaE = def_quantaE;            
            density = def_density;
            xRes    = def_xRes;
            xLims   = def_xLims;      
            optFilter = def_optFilter;
        elseif nargin == 3
            quantaE = def_quantaE;
            density = def_density;
            xRes    = def_xRes;
            xLims   = def_xLims;      
            optFilter = def_optFilter;
        elseif nargin == 4
            density = def_density;
            xRes    = def_xRes;
            xLims   = def_xLims;
            optFilter = def_optFilter;
        elseif nargin == 5
            xRes    = def_xRes;
            xLims   = def_xLims;
            optFilter = def_optFilter;
        elseif nargin == 6
            xLims   = def_xLims;
            optFilter = def_optFilter;
        elseif nargin == 7
            optFilter = def_optFilter;
        elseif nargin == 8
            % nothing to be done
        else
            errordlg('Too many inputs')
        end

    %% Main IMPLEMENTATION

        % create the wavelength vector
            lambda = (xLims(1) : xRes : xLims(2)); % [nm]      

            % whos 
            
        % Constants        

            % notice that the majorcase / minorcase notations are
            % reversed in the paper of Govardovskii compared to the
            % original Lamb paper!

            % modified the most from Lamb
            A = 0.8795 + (0.0459 * exp(-1 * (((peak_nm - 300)^2) / 11940))); % Lamb: 0.880

                % for peak_nm > 500 nm, modified equation keeps a
                % virtually constant, but as peak_nm decreases below
                % that, a starts growing progressively, increasing the
                % long-wave slope

            % only slight modifications of the constants to acquire
            % better fit
            B = 0.922;  % Lamb: 0.924
            C = 1.104;  % Lamb: 1.104
            D = 0.674;  % Lamb: 0.655

            a = 69.7;   % Lamb: 70
            b = 28.0;   % Lamb: 28.5
            c = -14.9;  % Lamb: -14.1

        % Vectorize the variables to avoid using a for-loop
            A = A * ones(1, length(lambda));
            B = B * ones(1, length(lambda));
            C = C * ones(1, length(lambda));
            D = D * ones(1, length(lambda));

            a = a * ones(1, length(lambda));
            b = b * ones(1, length(lambda));
            c = c * ones(1, length(lambda));

        % Defining denominator terms as their own variables to make the
        % code a bit more readable              
            denom_term1 = exp(a .* (A - (peak_nm ./ lambda) ) );
            denom_term2 = exp(b .* (B - (peak_nm ./ lambda) ) );
            denom_term3 = exp(c .* (C - (peak_nm ./ lambda) ) );
            denom_term4 = D;         

        %% Calculating spectral sensitivities 
        
            % (i.e. Lamb nomogram, ALPHA band)
        
                % Equation 2' from Lamb (1995)    
                S_alpha_initLinQ = (1 ./ (denom_term1 + denom_term2 + denom_term3 + denom_term4))';

            % BETA band
            
                peakBeta = 189 + (0.315 * peak_nm); % nm
                b        = -40.5 + (0.195 * peak_nm); % nm
                A_beta   = 0.26;
                x        = lambda - peakBeta;

                % the beta band absorption template
                S_betaBand_initLinQ = (A_beta * exp(-1* ((x/b) .^2)))';
            
        %% DETERMINE THE OUTPUT depending on the input FLAGs given
        S =  nomog_defineOutputCommon(S_alpha_initLinQ, S_betaBand_initLinQ,...
             bands, linLog, quantaE, density, lambda', optFilter);
        
        