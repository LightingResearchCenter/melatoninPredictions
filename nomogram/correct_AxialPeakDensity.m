function sensitivity_LOGcorr = correct_AxialPeakDensity(sensitivityIn, density)

     % for the correction equation, see Lamb (1995) for example            
     % Lamb (1995) used 0.27 for cones and 0.40 for rods

     % input in LOG units
     % ==================

     % for quick testing / debugging
     if nargin == 0
        load nomogramTest.mat
        density = 1;
        l = (380:0.1:780)';
        close all        
     else
        %l = (380:0.1:780)';
     end

     
     % create a vector of ones to match the length of the spectral
     % sensitivity
     one = ones(length(sensitivityIn),1); % unit vector

     % corrected, in LINEAR units
     transmittance = (10 .^ (-1 * (density * (10 .^ sensitivityIn))));     
     sensitivity_LINcorr = one - transmittance;
     
     % corrected, in LOG units
     sensitivity_LOGcorr = log10(sensitivity_LINcorr);     
     
     % debug
     %{
     figure
         subplot(2,1,1)
         plot(l,sensitivityIn, 'r', l,sensitivity_LOGcorr, 'b')
         legend('In LOG', 'Out LOG')
         
         subplot(2,1,2)
         plot(l,sensitivity_LINcorr, 'g', l, transmittance, 'r')
         legend('Out LIN ', 'Pigment Transmittance')
     %}
     