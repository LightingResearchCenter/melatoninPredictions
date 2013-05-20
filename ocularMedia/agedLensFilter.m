function lensFilter = agedLensFilter(age, lambda, offset)

    % WRAPPER Function to create a LENS density template as a function of
    % input, now only uses the van de Kraats and van Norren (2007) model   
    
    % INPUT
    %   age        - scalar defining the "age" of the ocular media
    %   lambda     - wavelength vector for the desired filter (e.g. 380:1:780)
    %   offset     - neutral offset to longer wavelengths (i.e. how much is
    %                light filtered in the long wavelengths where there is no
    %                age-dependency (0.111 in original paper)
    
    % OUTPUT
    %   lensFilter - spectral tranmittance of the LENS, LINEAR
    
    % call the subfunction
    templates.lens.vanDeKraats2007 = lensModel_vanDeKraats2007(age, lambda, offset);                
    
    % assign the sum of all components as the simple filter
    lensFilter = templates.lens.vanDeKraats2007.totalMedia;
    
    % convert from LOG DENSITY to LINEAR TRANSMITTANCE
    lensFilter = 10 .^ (lensFilter * -1);

function lensMedia = lensModel_vanDeKraats2007(A, lambda, offset)

    % if lambda is saved column-wise then transpose the input
    [m,n] = size(lambda);
    if n > m
        lambda = lambda';
    end

    % van de Kraats and van Norren (2004) - TRYPTOPHAN
    ones273 = 273 * ones(length(lambda),1);    
    lensMedia.TP = 14.19 * 10.68 * exp(-1 .* ((0.057 * (lambda - ones273))) .^ 2);

    % van de Kraats and van Norren (2004) - YOUNG LENS
    ones370 = 370 * ones(length(lambda),1);
    lensMedia.LY = (0.998 - (0.000063 * (A .^ 2))) * 2.13 * exp(-1 .* ((0.029 * (lambda - ones370))) .^ 2);

    % van de Kraats and van Norren (2004) - OLD LENS UV
    ones325 = 325 * ones(length(lambda),1);
    lensMedia.LOUV = (0.059 + (0.000186 * (A .^ 2))) * 11.95 * exp(-1 .* ((0.021 * (lambda - ones325))) .^ 2);

    % van de Kraats and van Norren (2004) - OLD LENS
    lensMedia.LO = (0.016 + (0.000132 * (A .^ 2))) * 1.43 * exp(-1 .* ((0.008 * (lambda - ones325))) .^ 2);                   

    % assign the scatter again
    scatterStruct = calc_rayleighScatter(A, lambda);
    lensMedia.rayleighScatter = scatterStruct.rayleighKraatsNorren2007_largeFields;

    % van de Kraats and van Norren (2004) - TOTAL MEDIA,
    % lenses + RAYLEIGH SCATTER + neutral offset [Equation
    % (8) of the paper]                    
    lensMedia.totalMedia = lensMedia.rayleighScatter...
                            + lensMedia.TP...
                            + lensMedia.LY...
                            + lensMedia.LOUV...
                            + lensMedia.LO...
                            + offset;

function scatterStruct = calc_rayleighScatter(A, lambda)
            
    % van de Kraats and van Norren (2007) - RAYLEIGH SCATTER
    % http://dx.doi.org/10.1364/JOSAA.24.001842
    % + http://dx.doi.org/10.1016/j.exer.2005.09.007
    
    offs = 0.446;
    scatterStruct.rayleighKraatsNorren2007_1deg = (offs + (0.000031 * (A .^ 2))) * ((400 ./ lambda) .^ 4);
    
    offs = 0.225;
    scatterStruct.rayleighKraatsNorren2007_largeFields = (offs + (0.000031 * (A .^ 2))) * ((400 ./ lambda) .^ 4);

    % dlmwrite(['rayleighScatter_age', num2str(handles.age), '_vanDeKraats_380_to_780nm_LIN.txt'],...
    %           [lambda 10 .^ (-1 * scatterStruct.rayleighKraatsNorren2007_largeFields)], 'Delimiter', '\t')



    
