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
