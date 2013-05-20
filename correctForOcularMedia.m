function [filt, lensFilter] = correctForOcularMedia(lambdaSimSPD, SPD, irrad, age, lambda, offset, path); 

    % NOTE!
    % Now all lensFilters refer to the relative transmittance
    % compared to the 25 yr old observer, so it is possible to
    % have transmittance values that are larger than 1

        %% Aphakic
        j = 1;
        cd(path.ocularMedia)
        lensFilterStd = agedLensFilter(age(j), lambda, offset); % 25 yr std observer
        cd(path.mainCode)

            % no the aphakic observer would have an increase of
            % retinal irradiance, so if we divide the original
            % light SPD with the standard observer we get the
            % retinal irradiance that would occur without the human
            % crystalline lens
            % go through the lights
            lensFilter{j} = 1 ./ lensFilterStd; % larger than 1 (transmittance)

            for i = 1 : length(lambdaSimSPD) 
                for k = 1 : length(irrad)
                    filt{j}{i}{k} = SPD.orig{1}{i}{k} .* lensFilter{j};
                end
            end


        %% 25yr old
        j = 2;
        % The same as the original model, no computations
        % needed

            cd(path.ocularMedia)
            lensFilter{j} = agedLensFilter(age(j-1), lambda, offset); % 25 yr std observer
            lensFilter{j} = lensFilter{j} ./ lensFilterStd; % normalize
                % the original model was for 25yr standard observer
                % so no additional lens filtering is needed, and
                % the phakic and 65 yr. old observers are related
                % to the zero attenuation
            cd(path.mainCode);

            for i = 1 : length(lambdaSimSPD) 
                for k = 1 : length(irrad)
                    filt{j}{i}{k} = SPD.orig{1}{i}{k} .* lensFilter{j};
                end
            end

        % 65yr old
        j = 3;
        cd(path.ocularMedia)
        lensFilter{j} = agedLensFilter(age(j-1), lambda, offset);
        cd(path.mainCode)

            % no the 65 yr old observer would have attenuated
            % retinal irradiance by the amount specified in
            % lensFilterTmp compared to aphakic observer, but as we
            % want the increased attenuation compared to the 25 yr
            % old observer
            lensFilter{j} = lensFilter{j} ./ lensFilterStd; % smaller than 1 (transmittance)

            for i = 1 : length(lambdaSimSPD) 
                for k = 1 : length(irrad)
                    filt{j}{i}{k} = SPD.orig{1}{i}{k} .* lensFilter{j};    
                end
            end 

        %                 figure
        %                 subplot(1,2,1)
        %                 plot(lambda, lensFilter{1})
        %                 xlim([400 600])
        %                 subplot(1,2,2)
        %                 plot(lambda, lensFilter{2}, 'r', lambda, lensFilter{3}, 'b')
        %                 xlim([400 600])
        %                 pause  