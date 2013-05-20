function melanopic = create_MelanopicCurve(lambda, lensFilter, templatePath)

    % import the data
    dataRaw = importdata(fullfile(templatePath, 'melanopic_humans_380to780_1nm_steps.txt'));
        wavelength = dataRaw.data(:,1);
        response = dataRaw.data(:,2);
        
    % weigh with the lens filter
    melanopic = response .* lensFilter;
    
    