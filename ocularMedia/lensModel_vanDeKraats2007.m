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