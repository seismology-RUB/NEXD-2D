% Script to print output from calc_vel.m
% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

if number_of_fluids == 2
    if ismember(used_model,[1,3,5]) == true
        fprintf('Used permeability model: van Genuchten (1980) \n');
    elseif ismember(used_model,[2,4,6]) == true
        fprintf('Used permeability model: Brooks & Corey (1964) \n');
    end
end
if number_of_fluids == 0 || phi_in == 0
    fprintf('\n');
    fprintf('Wavespeed of P for material %g: %g m/s \n', u, max(v1));
    fprintf('Wavespeed of S for material %g: %g m/s \n', u, max(v2));
elseif number_of_fluids > 0
    fprintf('\n');
    fprintf('--------------------------------------------------- \n');
    fprintf('Maximal wavespeed of P1 for material %g: %g m/s \n', u, max(v1));
    fprintf('Minimal wavespeed of P1 for material %g: %g m/s \n', u, min(v1));
    fprintf('Wavespeed at the most energetic frequency of P1: %g m/s \n', v1_fc(2));
    fprintf('--------------------------------------------------- \n');
    fprintf('Maximal wavespeed of S for material %g: %g m/s \n', u, max(v2));
    fprintf('Minimal wavespeed of S for material %g: %g m/s \n', u, min(v2));
    fprintf('Wavespeed at the most energetic frequency of S: %g m/s \n', v2_fc(2));
    if attenuation == true
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal inverse quality factor of P1 for material %g: %.15g \n', u, max(Q1));
        fprintf('Minimal inverse quality factor of P1 for material %g: %.15g \n', u, min(Q1));
        fprintf('Inverse quality factor at the most energetic frequency of P1: %.15g \n', Q1_fc(2));
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal inverse quality factor of S for material %g: %.15g \n', u, max(Q2));
        fprintf('Minimal inverse quality factor of S for material %g: %.15g \n', u, min(Q2));
        fprintf('Inverse quality factor at the most energetic frequency of S: %g \n', Q2_fc(2));
    end
    if number_of_fluids == 1
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal wavespeed of P2 for material %g: %g m/s \n', u, max(v3));
        fprintf('Minimal wavespeed of P2 for material %g: %g m/s \n', u, min(v3));
        fprintf('Wavespeed at the most energetic frequency of P2: %g m/s \n', v3_fc(2));
        if attenuation == true
            fprintf('--------------------------------------------------- \n');
            fprintf('Maximal inverse quality factor of P2 for material %g: %.15f \n', u, max(Q3));
            fprintf('Minimal inverse quality factor of P2 for material %g: %.15f \n', u, min(Q3));
            fprintf('Inverse quality factor at the most energetic frequency of P2: %.15f \n', Q3_fc(2));
        end
    elseif number_of_fluids == 2
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal wavespeed of P2 for material %g: %g m/s \n', u, max(v3));
        fprintf('Minimal wavespeed of P2 for material %g: %g m/s \n', u, min(v3));
        fprintf('Wavespeed at the most energetic frequency of P2: %g m/s \n', v3_fc(2));
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal wavespeed of P3 for material %g: %g m/s \n', u, max(v4));
        fprintf('Minimal wavespeed of P3 for material %g: %g m/s \n', u, min(v4));
        fprintf('Wavespeed at the most energetic frequency of P3: %g m/s \n', v4_fc(2));
        if attenuation == true
            fprintf('--------------------------------------------------- \n');
            fprintf('Maximal inverse quality of P2 for material %g: %.15f \n', u, max(Q3));
            fprintf('Minimal inverse quality of P2 for material %g: %.15f \n', u, min(Q3));
            fprintf('Inverse quality at the most energetic frequency of P2: %.15f \n', Q3_fc(2));
            fprintf('--------------------------------------------------- \n');
            fprintf('Maximal inverse quality of P3 for material %g: %.15f \n', u, max(Q4));
            fprintf('Minimal inverse quality of P3 for material %g: %.15f \n', u, min(Q4));
            fprintf('Inverse quality at the most energetic frequency of P3: %.15f \n', Q4_fc(2));
        end
    end
end
if number_of_fluids >= 1 
    fprintf(' \n');
    fprintf('Maximal wavelength of P1 for material %g: %.6f m \n', u, max(wavelength1));
    fprintf('Minimal wavelength of P1 for material %g: %.6f m \n', u, min(wavelength1));
    fprintf('--------------------------------------------------- \n');
    fprintf('Maximal wavelength of S for material %g: %.6f m \n', u, max(wavelength2));
    fprintf('Minimal wavelength of S for material %g: %.6f m \n', u, min(wavelength2));
    fprintf('--------------------------------------------------- \n');
    fprintf('Maximal wavelength of P2 for material %g: %.6f m \n', u, max(wavelength3));
    fprintf('Minimal wavelength of P2 for material %g: %.6f m \n', u, min(wavelength3));
    if number_of_fluids == 2
        fprintf('--------------------------------------------------- \n');
        fprintf('Maximal wavelength of P3 for material %g: %.6f m \n', u, max(wavelength4));
        fprintf('Minimal wavelength of P3 for material %g: %.6f m \n', u, min(wavelength4));
    end
end
