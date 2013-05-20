function S = nomog_defineOutputCommon(alpha, beta, bands, linLog, quantaE, density, lambda, optFilter)   

    %% quick fix for the common path
    currDir = pwd;
    cd ..;
    libRoot     = pwd;    
    commonPath  = fullfile(libRoot, 'common');
    cd(currDir)

    %% the wanted bands
    if strcmp(bands, 'both') == 1
        S = alpha + beta;
    elseif strcmp(bands, 'alpha') == 1
        S = alpha;
    elseif strcmp(bands, 'beta') == 1
        S = beta;
    end    
    
    %% whether in QUANTA or in ENERGY
    if strcmp(quantaE, 'Q') == 1
        % no changes needed
    elseif strcmp(quantaE, 'E') == 1       
        cd(commonPath)
        S = convert_fromQuantaToEnergy(S, lambda);
        S = S / max(S);
        cd(currDir)
    end
    
    %% whether in LIN or in LOG
    if strcmp(linLog, 'lin') == 1
        % no changes needed
        
    elseif strcmp(linLog, 'log') == 1
        S = log10(S);
    end  
    
    %% whether nomogram is corrected for axial density
    if density ~= 0        
        
        % convert to LOG
        if strcmp(linLog, 'lin') == 1
            S = log10(S);
        end
        
        % now the input S is in LOG units, as well as the output        
        S = correct_AxialPeakDensity(S, density);
        
        % normalize or not, no boolean switch now defined, so we always
        % normalize
        % S = S / max(S);
        
        % put back to LINEAR domain
        if strcmp(linLog, 'lin') == 1
            S = 10 .^ S;
        end
        
    else

    end
    
    
    %% if the output want to be filtered with LENS, MACULAR PIGMENT or 
    % any kind of optical filter, the optical filter is given as a template
    % in DENSITY and in LOG units
    if ~isempty(optFilter)
       
        if strcmp(linLog, 'lin') == 1
            S = log10(S); % convert to LOG
            S = S - optFilter;
            S = 10 .^ S; % convert back to LIN
        elseif strcmp(linLog, 'log') == 1
            S = S - optFilter;
        end    
        
    end
