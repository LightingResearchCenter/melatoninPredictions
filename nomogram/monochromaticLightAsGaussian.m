function SPD = monochromaticLightAsGaussian(peak, FWHM, xRes, xLimits)

    % to demo how to reconstruct the used light in some experiment based on
    % the peak wavelength and the half-bandwidth (HBW or FWHM)
    % Petteri Teikari, INSERM, 2010, petteri.teikari@inserm.fr
    
    %% user-defined parameters
    %peak = 480; % nm
    %FWHM = 10; % +/- nm 
    
        
    
        % sigma / deviation of the gaussian is defined as a function of FWHM
        sigma = (FWHM*2) / 2.355; % e.g., http://cp.literature.agilent.com/litweb/pdf/5980-0746E.pdf
    
        % create a wavelength vector        
        lambda = (xLimits(1) : xRes : xLimits(2))';
    
    %% create a GAUSSIAN Spectral Power Distribution (SPD) with the above
    %% parameters and lambda
    f = gauss_distribution(lambda, peak, sigma);
    SPD = f / max(f); % normalize       
    

        
    % Get the HBW / FWHM of the created light source, just for debug
    % hbw = getHBW(lambda, f);
   
    
    function f = gauss_distribution(x, mu, s)
    %% subfunction for Gaussian distribution
    
        % x         x vector (i.e. wavelength)
        % mu        mean / peak nm
        % sigma     standard deviation
        
        p1 = -.5 * ((x - mu)/s) .^ 2;
        p2 = (s * sqrt(2*pi));
        f = exp(p1) ./ p2;  
    
    
    function hbw = getHBW(x, y)  
    %% Subfunction to calculate the HBW of the created light
        
        % output
            % hbw           structure containing various parameters
        % inputs
            % x             wavelength of the light SPD
            % y             spectral irradiance/radiance/photon flux/etc. of the light
        
        hbw_safeLimit = 0.03; % safe limit to find the peak
        hbw_upsampleFactor = 10; % upsample factor for better spectral res.
        hbw_interpMode = 'pchip'; % interpolation method for upsampling

        % get the half-bandwidth of the spectrum (assuming that it
        % is a monochromatic spectrum with a Gaussian form)

            % normalize the spectrum
            y = y / max(y);

            % get a small portion around the peak for further interpolation 
            ind_min = find(y > hbw_safeLimit, 1, 'first');
            ind_max = find(y > hbw_safeLimit, 1, 'last');

                x_hbw = x(ind_min:ind_max);
                y_hbw = y(ind_min:ind_max);                

                vectorSizeReduction = (length(x_hbw) / length(x)); % debug variable                    

            % interpolate the vectors for better spatial resolution
            noOfPoints = (length(x_hbw) * hbw_upsampleFactor);
            x_interp = linspace(min(x_hbw), max(x_hbw), noOfPoints)';
            pp = interp1(x_hbw, y_hbw, hbw_interpMode, 'pp');                
            y_interp = ppval(pp, x_interp);               

            % get the difference between the interpolated y and 0.5
            % (brute force method) to find the points closest to 0.5
            % (in other words for the hbw points), it is very likely
            % that the y vector does not contain values EXACTLY 0.5000
            y_diff = abs(y_interp - 0.5);                

            % sort the rows in the hope that first two values would
            % correspond to the HBW points, this could be improved at
            % some point adding more constraints if it seems that it
            % starts to fail
            [y_diff,index] = sortrows(y_diff);      

                index_1 = index(1); index_2 = index(2);
                wavel_1 = x_interp(index_1); wavel_2 = x_interp(index_2);

                hbw.min = min([wavel_1 wavel_2]);
                hbw.max = max([wavel_1 wavel_2]);                                     
                [maxValue maxInd] = max(y_interp); 
                hbw.peak = x_interp(maxInd);
                hbw.range = hbw.max - hbw.min;

                hbw.interp_x = x_interp;
                hbw.interp_y = y_interp;
                hbw.wavel_1 = wavel_1;
                hbw.wavel_2 = wavel_2;
                hbw.index1 = index_1;
                hbw.index2 = index_2;

                hbw.peakRounded = round(hbw.peak);