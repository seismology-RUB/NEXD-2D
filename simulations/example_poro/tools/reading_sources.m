% This script calculates the highest frequency for each wavalet type
% with a fourier transformation.
%
% The following wavelet types are used:
% Green:    1
% Ricker:   2
% Sin:      3
% External: 4
%
% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

cd (workpath)
cd ../data
figure_number = figure_number+src;
f0            = frequency_vec(src); % Highest and most energetic frequency (in Hertz) for the Ricker wavelet
omega_p       = 2*pi*f0;            % Most energetic angular frequency (in radians per second)
wavelet_type  = stf(src);


if wavelet_type ~= 1 && wavelet_type ~= 2 && wavelet_type ~= 3 && wavelet_type ~= 4
    fprintf('Unknown type of source time function. \n');
    cd (workpath)
    cd ../tools
    return
elseif wavelet_type == 1
    omega_green = linspace(0,25*f0,1e6);
    Fourier_green = (exp(-omega_green.^2 ./ (4*pi^2*f0^2))) / (sqrt(2*pi)); % analytical fourier transformation of the Green wavelet

    % Finding roots of the green fourier transform
    roots_Fourier_green = find(Fourier_green <= wavelet_amplitude*max(Fourier_green));
    omega_max = omega_green(roots_Fourier_green(1));

elseif wavelet_type == 2
    omega_ricker   = linspace(0,25*f0,1e6);
    Fourier_ricker = (2.*omega_ricker.^2) ./ (sqrt(pi)*omega_p^3) .* (exp(-(omega_ricker.^2) ./(omega_p^2))); % analytical fourier transformation of the Ricker wavelet

    % Finding roots of F for maximal value of the angular frequency
    roots         = find(Fourier_ricker <= wavelet_amplitude*max(Fourier_ricker));
    droots        = diff(roots);
    discontinuty  = find(droots==max(droots))+1;
    discontinuty1 = roots(discontinuty);
    omega_max     = omega_ricker(discontinuty1);

elseif wavelet_type == 3
    omega_max = 2*pi*f0;

elseif wavelet_type == 4
    external_table = str2num(char(table2array(readtable(char(extwavelet{src})))));  % Reading the .txt file which includes the external wavelet
    time           = external_table(:,1);                                           % Time vector of the external wavelet
    amplitude      = external_table(:,2);                                           % Amplitude vector of the external wavelet

    sampling_rate  = length(amplitude)/(max(time)-min(time));                       % sampling rate calculated by the number of samples and the length of the external wavelet

    trans_length  = length(amplitude);                                              % Transformation of the length
    fast_ft       = fft(amplitude,trans_length);                                    % Fast Fourier Transform
    frequency     = (0:trans_length-1)*(sampling_rate/trans_length);                % Frequency range
    power         = fast_ft.*conj(fast_ft)/trans_length;                            % Power of the fourier transform

    % Parameters for roots of external wavelet
    angular_frequency_roots = 2*pi.*frequency(1:floor(trans_length/2));             % Frequency in rad/s
    power_roots             = power(1:floor(trans_length/2));                       % power of the external wavelet

    % Finding roots of fourier transformation for maximal value of the angular frequency
    roots_Fourier_external                      = find(power_roots <= wavelet_amplitude*max(power_roots));
    derivative_of_roots_Fourier_external        = diff(roots_Fourier_external);
    if isempty(find(derivative_of_roots_Fourier_external>1, 1))
        find_discontinuty = 1;
    else
        find_discontinuty                       = find(derivative_of_roots_Fourier_external==max(derivative_of_roots_Fourier_external))+1;
    end
    find_entry_in_roots_Fourier_external        = roots_Fourier_external(find_discontinuty);
    omega_max                                   = angular_frequency_roots(find_entry_in_roots_Fourier_external); % Determination of maximal value of the angular frequency
end

% Printing what wavelet type is used and the center frequency f0
fprintf('\n');
if wavelet_type == 1
    fprintf('Source-time-function at source %g: Green \n',src);
elseif wavelet_type == 2
    fprintf('Source-time-function at source %g: Ricker \n',src);
elseif wavelet_type == 3
    fprintf('Source-time-function at source %g: sin^3 \n',src);
elseif wavelet_type == 4
    fprintf('Source-time-function at source %g: External \n',src);
end
if wavelet_type ~= 4
    fprintf('Center frequency of source-time-function: %g Hz \n',f0);
end
