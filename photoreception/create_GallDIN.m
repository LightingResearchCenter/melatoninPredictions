function gallDIN = create_GallDIN(lambda, lensFilter, templatePath)

    % import the data
    dataRaw = importdata(fullfile(templatePath, 'cLambda-Philips_LIN_380to780nm_1nm.txt'));
        wavelength = dataRaw.data(:,1);
        response = dataRaw.data(:,2);
        
    % weigh with the lens filter
    gallDIN = response .* lensFilter;
    
    