% Imports the Corneal sensitivities of different photoreceptors from disk
function S = import_CornealSensitivities(path)

    cd(path)

    % CONE fundamentals
    fileName = 'coneFundamentals_LINss10deg_380to780nm_1nm_Lin-E.txt';    
    dataRawCones = importdata(fileName, '\t', 1);

    % V(lambda), photopic
    fileName = 'v_lambda_linCIE2008v10e_LIN_380to780nm_1nm.txt';    
    dataRawPhotop = importdata(fileName, '\t', 1);

    % V'(lambda), scotopic
    fileName = 'v_lambda_Scotopic_1951e_LIN_380to780nm_1nm.txt';    
    dataRawScot = importdata(fileName, '\t', 1);

    % Melanopic, (Lucas et al. 2011)
    fileName = 'melanopic_humans_380to780_1nm_steps.txt';    
    dataRawMelanop = importdata(fileName, '\t', 1);
    
    % separate to variables   
    
        % Lambda vector
        S.lambda = dataRawCones.data(:,1);

        % SWS
        S.SWS = dataRawCones.data(:,4);

        % MWS
        S.MWS = dataRawCones.data(:,3);

        % LWS
        S.LWS = dataRawCones.data(:,2);

        % V(l)
        S.Vl = dataRawPhotop.data(:,2);

        % V'(l)
        S.VlPrime = dataRawScot.data(:,2);

        % Melanopic
        S.Melanop = dataRawMelanop.data(:,2);

        % get peaks
        names = fieldnames(S);
        peakTemplate = zeros((length(names)-1),1);
        for k = 2 : length(names) % skip wavelength       
            [~,ind] = max(S.(names{k}));
            peakTemplate(k-1) = S.lambda(ind);                
        end

        S.headers = {['SWS (', num2str(peakTemplate(1)), ' nm)'];...
                      ['MWS (', num2str(peakTemplate(2)), ' nm)'];...
                      ['LWS (', num2str(peakTemplate(3)), ' nm)'];...
                      ['V(\lambda) (', num2str(peakTemplate(4)), ' nm)'];...
                      ['V''(\lambda) (', num2str(peakTemplate(5)), ' nm)'];...
                      ['Melanopic (', num2str(peakTemplate(6)), ' nm)']};