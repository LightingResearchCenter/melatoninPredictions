%Brainard and Thappan Data
% [wavelength  response]
BandT = [420 0.256; ...
         424 0.815; ...
         440 0.953; ...
         456 1.000; ...
         460 1.000; ...
         472 0.8560; ...
         480 0.7259; ...
         496 0.5869; ...
         505 0.7916; ...
         520 0.5202; ...
         530 0.3958; ...
         548 0.1428; ...
         555 0.1089; ...
         575 0.0554; ...
         600 0.0282];
waveBT = BandT(:,1);
BT = BandT(:,2);

wave = (380:1:780)';
criterion = 300;
ICLA = zeros(length(wave),1);
for j = 1:length(wave)
    spd = zeros(length(wave),1);
    % Uncomment lines 24 to 29 (or51 for plots) and comment out line 31 for Gaussian shaped quasi-monochromatic stimuli 
    %FWHM = 15; %Full width half maximum (nm)
    %spd = exp(-1/(2*FWHM)*(wave-(wave(j))).^2);
    %spd = spd/sum(spd); % normalize to a radiant power of 1
    %figure(1)
    %plot(wave,spd)
    %hold on
    
    spd(j) = 1;
    steps = -4:0.2:3; % steps of log10(irradiance) (W/m^2)
    CLA = zeros(length(steps),1);
    I = zeros(length(steps),1);
    for i = 1:length(steps)
        I(i) = 10^steps(i); % irradiance (W/m^2)
        P = spd*I(i); % scaled spd
 
        CLA(i) = CLA_postBerlinCorrMelanopsin_02Oct2012([wave P]);
        
    end
    CLA
    ICLA(j) = interp1(CLA,I,criterion,'linear'); % find the irradiance that gives the criterion response
    display(['Wavelength = ' num2str(wave(j),'%d')]);
end
ICLAeff = 1./ICLA; % efficiency  = 1/(criterion response)
ICLAeff = ICLAeff/max(ICLAeff); % normalize to a maximum of 1

%figure(1)
%hold off
figure(2)
plot(wave,ICLAeff,'k--')
axis([400 650 0 1.2])
hold on
plot(waveBT,BT,'rd')
legend('Model','Brainard and Thapan')

save('ReaCurve.mat', 'wave', 'ICLAeff')

