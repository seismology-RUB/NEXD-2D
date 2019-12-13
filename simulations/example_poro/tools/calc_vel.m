% Tool to calculate wavespeeds for porousmaterials in two dimensions
% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

clear all;
warning('off','all')
workpath = pwd;
cd (workpath)
cd ../
mkdir out
cd out/

%%###############################
%% CONVENIENCE PARAMETERS
plot_diag           = true;       % if plot_diag = true, the plots are displayed
save_plot           = false;      % if saves_plot = true, the plot is saved in available formattypes in matlab
formattype          = '-dpng';
attenuation_default = true;       % if attenuation_default = true, the invers quality factor is calculated
logfile             = true;       % Set true if logfile should be created, otherwise false
k_default           = 40;         % Length of all vectors
wavelet_amplitude   = 0.001;      % minimal value of amplitude of fourier transformation of the wavelet,
                                  % used for frequency estimation. E.g. wavelet_amplitude = 0.001, then 0.1 % of the
                                  % highest amplitude for the highest frequency
%%###############################
%% Logfile
if logfile
    diary on
    delete('logfile_calc_vel')
    diary('logfile_calc_vel')
end
save_plot_diag = plot_diag;
cd (workpath)
cd ../tools

%% READING INPUT FILES

run('reading_parfile')
sourcematrix = [frequency_vec, stf];
cmp_elements = true;

%% CALCULATE VELOCITIES FOR EACH SOURCE

figure_number = 0; % figure_number is important to create for each iteration a figure
for src = 1:nsrc
    cd (workpath)
    cd ../tools
    run('reading_sources')

    cd (workpath)
    cd ../tools
    run('proof_source')

    if cmp_elements == true | cmp_elements > 1

        % Allocate vectors for frequency stability
        s_frequ = zeros(1,number_of_materials);
        s1      = zeros(1,number_of_materials);
        s2      = zeros(1,number_of_materials);
        s3      = zeros(1,number_of_materials);
        s4      = zeros(1,number_of_materials);

        %% CALCULATE VELOCITIES FOR EACH MATERIAL
        for u = 1:number_of_materials
            number_of_fluids = number_of_fluids_save;
            cd (workpath)
            cd ../tools
            run('porousmaterial.m')

            if phi_in == 0
                number_of_fluids = 0;
            end

            % Calculation of further parameters for calcualting wavespeeds
            if number_of_fluids == 0 % 0 fluids
                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        if phi_in == 0
                            b     = zeros(1,k); % Biot coefficient
                            inv_N = zeros(1,k); % Biot modulus
                        elseif phi_in ~= 0
                            b     = 1 - K_d ./ K_s;
                            inv_N = b ./ K_s;
                        end
                    elseif ismember(used_model,[5,6]) == true
                        b = 1/2 - sqrt(1/4 - K_d.*inv_N);
                    end
                    M        = 1./inv_N;
                    lambda_u = K_d + M.*b.^2 - 2/3.*my; % first Lamé Parameter (undrained)
                end

                phi1          = zeros(1,k);
                phi2          = zeros(1,k);
                M             = zeros(1,k);
                M1            = M;
                M2            = M;
                Mtilde        = M;
                M1tilde       = M;
                M2tilde       = M;
                rho_star      = rhos;
                gamma1        = zeros(1,k);
                gamma2        = zeros(1,k);
                Lambda1_star  = zeros(1,k);
                Lambda2_star  = zeros(1,k);
                Lambda11tilde = zeros(1,k);
                Lambda21tilde = zeros(1,k);
                Lambda12tilde = zeros(1,k);
                Lambda22tilde = zeros(1,k);

            elseif number_of_fluids == 1 % 1 fluid
                if S1_in == 0
                    K1   = K2;
                    rho1 = rho2;
                    eta1 = eta2;
                end
                S1   = ones(1,k);
                S2   = zeros(1,k);
                phi1 = S1.*phi;
                phi2 = S2.*phi;
                rho  = (1-phi).*rhos + (phi1.*rho1 + phi2.*rho2);

                if ismember(used_model,[1,3,5]) == true   % van Genuchten
                    m = 1-1./fitting_n;
                    G = zeros(1,k);                       % dp_c/dS_1 van Genuchten
                elseif ismember(used_model,[2,4,6]) == true
                    G = zeros(1,k);                       % dp_c/dS_1 Brooks & Corey
                end

                K_f = (S1./K1).^(-1);                         % Bulk modulus fluid

                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        b     = 1 - K_d ./ K_s;                 % Biot coefficient
                        inv_N = (b-phi) ./ K_s;                 % inverse Biot modulus
                    elseif ismember(used_model,[5,6]) == true
                        b = (phi+1)/2 - sqrt(((phi+1).^2)/4 - phi - K_d.*inv_N);
                    end
                    M = (inv_N + phi./K_f).^(-1);
                    lambda_u = K_d + M.*b.^2 - 2/3.*my; % first Lamè Parameter (undrained)
                end

                omegac1 = eta1.*inv_T.*phi1 ./ (kappa.*rho1);
                omegac2 = eta2.*inv_T.*phi2 ./ (kappa.*rho2);
                kappa1  = kappa;
                kappa2  = kappa;

                % Implementing correction function
                C1 = ones(1,k);
                C2 = ones(1,k);
                if eta1_in ~= 0
                    C1 = sqrt(1 + 0.5i .* omega ./ omegac1);
                end
                if eta2_in ~= 0
                    C2 = sqrt(1 + 0.5i .* omega ./ omegac2);
                end

                if ismember(used_model,[1,2,5,6]) == true
                    M = (inv_N + phi./K_f).^(-1);
                end
                M1            = M.*(1+S1.*S2./K1 .* G);
                M2            = M;
                Mtilde        = M;
                M1tilde       = M.*((1-inv_N.*S1.^2./phi + S1./K1) .* G);
                M2tilde       = M;
                rho_star      = (1-phi).*rhos + (1-inv_T) .* (phi1.*rho1 + phi2.*rho2);
                gamma1        = phi1 .* inv_T;
                gamma2        = phi2 .* inv_T;
                Lambda1_star  = kappa1 ./ (C1.*phi1.*eta1);
                Lambda2_star  = zeros(1,k);
                Lambda11tilde = rho_star.*Lambda1_star ./ (phi1.*(inv_T.^2-inv_T));
                Lambda21tilde = zeros(1,k);
                Lambda12tilde = (inv_T./(rho1.*Lambda1_star) + 1./Lambda11tilde).^(-1);
                Lambda22tilde = zeros(1,k);

            elseif number_of_fluids == 2 % 2 fluids
                phi1 = S1.*phi;
                phi2 = S2.*phi;
                rho  = (1-phi).*rhos + (phi1.*rho1 + phi2.*rho2);

                % Calculate effective saturations
                S1eff = (S1-Sr1)./(1-Sr2-Sr1);
                S2eff = (S2-Sr2)./(1-Sr1-Sr2);

                if ismember(used_model,[1,3,5]) == true
                    m = 1-1./fitting_n;
                    %G = - rho1.*g./((fitting_n-1).*fitting_chi) .* (S1.^(-fitting_n./(fitting_n-1))-1).^(1./fitting_n-1).*S1.^(-(2.*fitting_n-1)./(fitting_n-1)); % dp_c(dS_1 van Genuchten
                    G = - (rho1.*g)./((fitting_n-1).*fitting_chi) .* ((S1eff.^(-1./m)-1).^(1./fitting_n-1)) .* S1eff.^(-(1./m+1)); % dp_c(dS_1 van Genuchten
                    %G = - (rho1.*g)./((fitting_n-1).*fitting_chi) .* ((S1.^(-1./m)-1).^(1./fitting_n-1)) .* S1.^(-(1./m+1)); % dp_c(dS_1 van Genuchten
                    %G = rho2.*g./((fitting_n-1).*fitting_chi) .* (1-S1).^(1-((2.*fitting_n-1)./(fitting_n-1))); % dp_c(dS_1 van Genuchten
                elseif ismember(used_model,[2,4,6]) == true
                    G = - p_b ./ (lambda_BC * (Sr2-Sr1)) .* S1eff.^(-(1+lambda_BC)./lambda_BC); % dp_c/dS_1 Brooks & Corey
                    %G = - p_b ./ lambda_BC .* S1.^(-(1+lambda_BC)./lambda_BC); % dp_c/dS_1 Brooks & Corey
                    %G = p_b ./ lambda_BC .* (1-S1).^(-(1+lambda_BC)./lambda_BC); % dp_c/dS_1 Brooks & Corey
                elseif ismember(used_model,[7,8,9]) == true
                    % Douglas Jr. et al.
                    G = - 2 .* A .* (1 ./ ((S1-Sr1).^3) + Sr2.^2 ./ ((1-S1).^3 .* (1-Sr1-Sr2).^2));
                end

                K_f = (S1./K1 + S2./K2).^(-1);  % Bulk modulus of both fluids

                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        b        = 1 - K_d ./ K_s;
                        inv_N    = (b-phi) ./ K_s;
                    elseif ismember(used_model,[5,6]) == true
                        b = (phi+1)/2 - sqrt(((phi+1).^2)/4 - phi - K_d.*inv_N);
                    end
                    M        = (inv_N + phi1./K1 + phi2./K2 - ((S1.*S2.*phi) ./ (K1.*K2) + S1.*S2.*inv_N .* (S1 ./ K2 + S2 ./ K1)).*G).^(-1);
                    lambda_u = K_d + (1+S1.*S2.*(S1./K2+S2./K1).*G).*M.*b.^2 - 2/3.*my;
                end

                % Calculation of realative permeability
                if ismember(used_model,[1,3,5]) == true
                    % van Genuchten
                    kappa1rel = sqrt(S1eff) .* (1 - (1-S1eff.^(1./m)).^m).^2;
                    kappa2rel = sqrt(1-S1eff) .* (1 - S1eff.^(1./m)).^(2.*m);
                elseif ismember(used_model,[2,4,6]) == true
                    % Brooks & Corey
                    kappa1rel = S1eff.^((2+3.*lambda_BC)./lambda_BC);
                    kappa2rel = (1-S1eff).^2 .* (1-S1eff.^((2+lambda_BC)./lambda_BC));
                elseif ismember(used_model,[7,8,9]) == true
                    % Douglas Jr. et al.
                    kappa1rel = (1 -  S1 ./(1-Sr2)).^2;
                    kappa2rel = ((S1-Sr1)./(1-Sr1)).^2;
                end

                omegac1 = eta1.*inv_T.*phi1 ./ (kappa.*rho1);
                omegac2 = eta2.*inv_T.*phi2 ./ (kappa.*rho2);
                kappa1  = kappa1rel .* kappa;
                kappa2  = kappa2rel .* kappa;

                % Implementing correction function
                C1 = ones(1,k);
                C2 = ones(1,k);
                if eta1_in ~= 0
                    C1 = sqrt(1 + 0.5i .* omega ./ omegac1);
                end
                if eta2_in ~= 0
                    C2 = sqrt(1 + 0.5i .* omega ./ omegac2);
                end

                M             = (inv_N + phi1./K1 + phi2./K2 - ((S1.*S2.*phi) ./ (K1.*K2) + S1.*S2.*inv_N .* (S1 ./ K2 + S2 ./ K1)).*G).^(-1);
                M1            = M .* (1 - S1.*S2./K1 .* G);
                M2            = M .* (1 - S1.*S2./K2 .* G);
                Mtilde        = M .* (1 + inv_N.*S1.*S2./phi .* G);
                M1tilde       = M .* (1 - (inv_N.*S1.^2./phi + S1./K1) .* G);
                M2tilde       = M .* (1 - (inv_N.*S2.^2./phi + S2./K2) .* G);
                rho_star      = (1-phi).*rhos + (1-inv_T) .* (phi1.*rho1 + phi2.*rho2);
                gamma1        = phi1 .* inv_T;
                gamma2        = phi2 .* inv_T;
                Lambda1_star  = kappa1 ./ (C1.*phi1.*eta1);
                Lambda2_star  = kappa2 ./ (C2.*phi2.*eta2);
                Lambda11tilde = (rho_star .* Lambda1_star) ./ (phi1.*(inv_T.^2-inv_T));
                Lambda21tilde = (rho_star .* Lambda2_star) ./ (phi2.*(inv_T.^2-inv_T));
                Lambda12tilde = (inv_T ./ (rho1.*Lambda1_star) + 1./Lambda11tilde).^(-1);
                Lambda22tilde = (inv_T ./ (rho2.*Lambda2_star) + 1./Lambda21tilde).^(-1);
            end

            % Set attenuation, depending on porosity and viscosity
            attenuation = attenuation_default;
            if number_of_fluids == 0
                attenuation = false;
            elseif number_of_fluids == 1 && (phi_in == 0 || norm(eta1) == 0)
                attenuation = false;
            elseif number_of_fluids == 2
                if phi_in == 0 || (norm(eta1) == 0 && norm(eta2) == 0 )
                    attenuation = false;
                end
            end

            %% Wavespeed calculation
            % Allocate wavespeed  and inverse quality factor vectors
            v1 = zeros(k,1);
            v2 = zeros(k,1);
            v3 = zeros(k,1);
            v4 = zeros(k,1);

            Q1 = zeros(k,1);
            Q2 = zeros(k,1);
            Q3 = zeros(k,1);
            Q4 = zeros(k,1);
            
            wavelength1 = zeros(k,1);
            wavelength2 = zeros(k,1);
            wavelength3 = zeros(k,1);
            wavelength4 = zeros(k,1);

            for  i = 1:k
                %% ELEMENTS OF MATRIX A
                A11   = 0;
                A12   = 0;
                A13   = 0;
                A14   = -(lambda_u(i)+2*my(i))+b(i)*(phi1(i)*M2(i)+phi2(i)*M1(i));
                A15   = 0;
                if number_of_fluids ~= 0
                    A16   = 0;
                    A17   = -phi1(i)*M2(i)*b(i);
                    A18   = 0;
                    A19   = 0;
                    A110  = -phi2(i)*M1(i)*b(i);
                    A1_11 = 0;
                end

                A21  = 0;
                A22  = 0;
                A23  = 0;
                A24  = -lambda_u(i)+b(i)*(phi1(i)*M2(i)+phi2(i)*M1(i));
                A25  = 0;
                if number_of_fluids ~= 0
                    A26  = 0;
                    A27  = -phi1(i)*M2(i)*b(i);
                    A28  = 0;
                    A29  = 0;
                    A210 = -phi2(i)*M1(i)*b(i);
                    A211 = 0;
                end

                A31  = 0;
                A32  = 0;
                A33  = 0;
                A34  = 0;
                A35  = -my(i);
                if number_of_fluids ~= 0
                    A36  = 0;
                    A37  = 0;
                    A38  = 0;
                    A39  = 0;
                    A310 = 0;
                    A311 = 0;
                end

                A41  = -1/rho_star(i);
                A42  = 0;
                A43  = 0;
                A44  = 0;
                A45  = 0;
                if number_of_fluids ~= 0
                    A46  = -gamma1(i)/rho_star(i);
                    A47  = 0;
                    A48  = 0;
                    A49  = -gamma2(i)/rho_star(i);
                    A410 = 0;
                    A411 = 0;
                end

                A51  = 0;
                A52  = 0;
                A53  = -1/rho_star(i);
                A54  = 0;
                A55  = 0;
                if number_of_fluids ~= 0
                    A56  = 0;
                    A57  = 0;
                    A58  = 0;
                    A59  = 0;
                    A510 = 0;
                    A511 = 0;
                end

                if number_of_fluids ~= 0
                    A61  = 0;
                    A62  = 0;
                    A63  = 0;
                    if phi_in ~= 0
                        A64  = M2(i)*b(i)-(phi1(i)*M2tilde(i)+phi2(i)*Mtilde(i));
                    elseif phi_in == 0
                        A64 = 0;
                    end
                    A65  = 0;
                    A66  = 0;
                    if phi_in ~= 0
                        A67  = phi1(i)*M2tilde(i);
                    elseif phi_in == 0
                        A67 = 0;
                    end
                    A68  = 0;
                    A69  = 0;
                    if phi_in ~= 0
                        A610 = phi2(i)*Mtilde(i);
                    elseif phi_in == 0
                        A610 = 0;
                    end
                    A611 = 0;

                    A71  = (inv_T(i) -1) / rho_star(i);
                    A72  = 0;
                    A73  = 0;
                    A74  = 0;
                    A75  = 0;
                    A76  = inv_T(i) / rho1(i) + (inv_T(i)*phi1(i)*(inv_T(i)-1)) / rho_star(i);
                    A77  = 0;
                    A78  = 0;
                    A79  = (inv_T(i)*phi2(i)*(inv_T(i)-1)) / rho_star(i);
                    A710 = 0;
                    A711 = 0;

                    A81  = 0;
                    A82  = 0;
                    A83  = (inv_T(i) -1) / rho_star(i);
                    A84  = 0;
                    A85  = 0;
                    A86  = 0;
                    A87  = 0;
                    A88  = 0;
                    A89  = 0;
                    A810 = 0;
                    A811 = 0;

                    A91  = 0;
                    A92  = 0;
                    A93  = 0;
                    A94  = M1(i)*b(i)-(phi1(i)*Mtilde(i)+phi2(i)*M1tilde(i));
                    A95  = 0;
                    A96  = 0;
                    if phi ~= 0
                        A97  = phi1(i)*Mtilde(i);
                    elseif phi_in == 0
                        A97 = 0;
                    end
                    A98  = 0;
                    A99  = 0;
                    A910 = phi2(i)*M1tilde(i);
                    A911 = 0;

                    A101  = (inv_T(i) -1) / rho_star(i);
                    A102  = 0;
                    A103  = 0;
                    A104  = 0;
                    A105  = 0;
                    A106  = (inv_T(i)*phi1(i)*(inv_T(i)-1)) / rho_star(i);
                    A107  = 0;
                    A108  = 0;
                    A109  = inv_T(i)/rho2(i) + (inv_T(i)*phi2(i)*(inv_T(i)-1)) / rho_star(i);
                    A1010 = 0;
                    A1011 = 0;

                    A11_1  = 0;
                    A112   = 0;
                    A113   = (inv_T(i) -1) / rho_star(i);
                    A114   = 0;
                    A115   = 0;
                    A116   = 0;
                    A117   = 0;
                    A118   = 0;
                    A119   = 0;
                    A1110  = 0;
                    A1111  = 0;
                end

                %% ELEMENTS OF MATRIX E
                E11   = 0;
                E12   = 0;
                E13   = 0;
                E14   = 0;
                E15   = 0;
                E16   = 0;
                E17   = 0;
                E18   = 0;
                E19   = 0;
                E110  = 0;
                E1_11 = 0;

                E21  = 0;
                E22  = 0;
                E23  = 0;
                E24  = 0;
                E25  = 0;
                E26  = 0;
                E27  = 0;
                E28  = 0;
                E29  = 0;
                E210 = 0;
                E211 = 0;

                E31  = 0;
                E32  = 0;
                E33  = 0;
                E34  = 0;
                E35  = 0;
                E36  = 0;
                E37  = 0;
                E38  = 0;
                E39  = 0;
                E310 = 0;
                E311 = 0;

                E41  = 0;
                E42  = 0;
                E43  = 0;
                if number_of_fluids == 0 || number_of_fluids == 2
                    E44  = -gamma1(i)/(rho_star(i)*Lambda1_star(i)) - gamma2(i)/(rho_star(i)*Lambda2_star(i));
                elseif number_of_fluids == 1
                    E44 = -gamma1(i)/(rho_star(i)*Lambda1_star(i));
                end
                E45  = 0;
                E46  = 0;
                E47  = gamma1(i)/(rho_star(i)*Lambda1_star(i));
                E48  = 0;
                E49  = 0;
                E410 = gamma2(i)/(rho_star(i)*Lambda2_star(i));
                E411 = 0;

                E51  = 0;
                E52  = 0;
                E53  = 0;
                E54  = 0;
                if number_of_fluids == 0 || number_of_fluids == 2
                    E55  = -gamma1(i)/(rho_star(i)*Lambda1_star(i)) - gamma2(i)/(rho_star(i)*Lambda2_star(i));
                elseif number_of_fluids == 1
                    E55 = -gamma1(i)/(rho_star(i)*Lambda1_star(i));
                end
                E56  = 0;
                E57  = 0;
                E58  = gamma1(i)/(rho_star(i)*Lambda1_star(i));
                E59  = 0;
                E510 = 0;
                E511 = gamma2(i)/(rho_star(i)*Lambda2_star(i));

                E61  = 0;
                E62  = 0;
                E63  = 0;
                E64  = 0;
                E65  = 0;
                E66  = 0;
                E67  = 0;
                E68  = 0;
                E69  = 0;
                E610 = 0;
                E611 = 0;

                E71  = 0;
                E72  = 0;
                E73  = 0;
                if number_of_fluids == 0 || number_of_fluids == 2
                    E74  = 1/Lambda12tilde(i) + 1/Lambda21tilde(i);
                elseif number_of_fluids == 1
                    E74 = 1/Lambda12tilde(i);
                end
                E75  = 0;
                E76  = 0;
                E77  = -1/Lambda12tilde(i);
                E78  = 0;
                E79  = 0;
                E710 = -1/Lambda21tilde(i);
                E711 = 0;

                E81  = 0;
                E82  = 0;
                E83  = 0;
                E84  = 0;
                if number_of_fluids == 0 || number_of_fluids == 2
                    E85  = 1/Lambda12tilde(i) + 1/Lambda21tilde(i);
                elseif number_of_fluids == 1
                    E85 = 1/Lambda12tilde(i);
                end
                E86  = 0;
                E87  = 0;
                E88  = -1/Lambda12tilde(i);
                E89  = 0;
                E810 = 0;
                E811 = -1/Lambda21tilde(i);

                E91  = 0;
                E92  = 0;
                E93  = 0;
                E94  = 0;
                E95  = 0;
                E96  = 0;
                E97  = 0;
                E98  = 0;
                E99  = 0;
                E910 = 0;
                E911 = 0;

                E101  = 0;
                E102  = 0;
                E103  = 0;
                E104  = 1/Lambda11tilde(i) + 1/Lambda22tilde(i);
                E105  = 0;
                E106  = 0;
                E107  = -1/Lambda11tilde(i);
                E108  = 0;
                E109  = 0;
                E1010 = -1/Lambda22tilde(i);
                E1011 = 0;

                E11_1  = 0;
                E112   = 0;
                E113   = 0;
                E114   = 0;
                E115   = 1/Lambda11tilde(i) + 1/Lambda22tilde(i);
                E116   = 0;
                E117   = 0;
                E118   = -1/Lambda11tilde(i);
                E119   = 0;
                E1110  = 0;
                E1111  = -1/Lambda22tilde(i);

                %% CREATING THE MATRICES
                if number_of_fluids == 0
                    A = [A11, A12, A13, A14, A15;
                         A21, A22, A23, A24, A25;
                         A31, A32, A33, A34, A35;
                         A41, A42, A43, A44, A45;
                         A51, A52, A53, A54, A55];

                    E = zeros(5);
                elseif number_of_fluids == 1
                    A = [A11, A12, A13, A14, A15, A16, A17, A18;
                         A21, A22, A23, A24, A25, A26, A27, A28;
                         A31, A32, A33, A34, A35, A36, A37, A38;
                         A41, A42, A43, A44, A45, A46, A47, A48;
                         A51, A52, A53, A54, A55, A56, A57, A58;
                         A61, A62, A63, A64, A65, A66, A67, A68;
                         A71, A72, A73, A74, A75, A76, A77, A78;
                         A81, A82, A83, A84, A85, A86, A87, A88];

                    E = [E11, E12, E13, E14, E15, E16, E17, E18;
                         E21, E22, E23, E24, E25, E26, E27, E28;
                         E31, E32, E33, E34, E35, E36, E37, E38;
                         E41, E42, E43, E44, E45, E46, E47, E48;
                         E51, E52, E53, E54, E55, E56, E57, E58;
                         E61, E62, E63, E64, E65, E66, E67, E68;
                         E71, E72, E73, E74, E75, E76, E77, E78;
                         E81, E82, E83, E84, E85, E86, E87, E88];

                elseif number_of_fluids == 2
                    A = [A11, A12, A13, A14, A15, A16, A17, A18, A19, A110, A1_11;
                         A21, A22, A23, A24, A25, A26, A27, A28, A29, A210, A211;
                         A31, A32, A33, A34, A35, A36, A37, A38, A39, A310, A311;
                         A41, A42, A43, A44, A45, A46, A47, A48, A49, A410, A411;
                         A51, A52, A53, A54, A55, A56, A57, A58, A59, A510, A511;
                         A61, A62, A63, A64, A65, A66, A67, A68, A69, A610, A611;
                         A71, A72, A73, A74, A75, A76, A77, A78, A79, A710, A711;
                         A81, A82, A83, A84, A85, A86, A87, A88, A89, A810, A811;
                         A91, A92, A93, A94, A95, A96, A97, A98, A99, A910, A911;
                         A101, A102, A103, A104, A105, A106, A107, A108, A109, A1010, A1011;
                         A11_1, A112, A113, A114, A115, A116, A117, A118, A119, A1110, A1111];

                     E = [E11, E12, E13, E14, E15, E16, E17, E18, E19, E110, E1_11;
                          E21, E22, E23, E24, E25, E26, E27, E28, E29, E210, E211;
                          E31, E32, E33, E34, E35, E36, E37, E38, E39, E310, E311;
                          E41, E42, E43, E44, E45, E46, E47, E48, E49, E410, E411;
                          E51, E52, E53, E54, E55, E56, E57, E58, E59, E510, E511;
                          E61, E62, E63, E64, E65, E66, E67, E68, E69, E610, E611;
                          E71, E72, E73, E74, E75, E76, E77, E78, E79, E710, E711;
                          E81, E82, E83, E84, E85, E86, E87, E88, E89, E810, E811;
                          E91, E92, E93, E94, E95, E96, E97, E98, E99, E910, E911;
                          E101, E102, E103, E104, E105, E106, E107, E108, E109, E1010, E1011;
                          E11_1, E112, E113, E114, E115, E116, E117, E118, E119, E1110, E1111];
                end

                %% WAVESPEED CALCUATION
                syms s % Creating symbolic variable for wavenumber

                if number_of_fluids == 0
                    order_of_identity_matrix = 5;
                elseif number_of_fluids == 1
                    order_of_identity_matrix = 8;
                elseif number_of_fluids == 2
                    order_of_identity_matrix = 11;
                end

                % Implementing fractional step
                if fractional_step == true && attenuation == false
                    E = zeros(order_of_identity_matrix);
                end

                if fractional_step == false
                    E1       = 1i*E;                                                    % Matrix E gets its imaginary part
                    equation = det(-omega(i)*eye(order_of_identity_matrix) + s*A+E1);   % Solution of the eigenvalue problem
                    solution = vpasolve(equation==0,s);                                 % Solves 'eqn' with respect to the wavenumber s
                elseif fractional_step == true
                    E1 = 1i*E;
                    equation     = det(-omega(i)*eye(order_of_identity_matrix) + s*A+E1);
                    solution_att = vpasolve(equation==0,s);

                    equation   = det(-omega(i)*eye(order_of_identity_matrix) + s*A);
                    solution_v = vpasolve(equation==0,s);
                end

                %% CALCUATION OF WAVESPEED AND INVERSE QUALITY FACTOR
                %% CREATING VECTORS WITH WAVESPEED AND INVERSE QUALITY FACTOR
                if fractional_step == false
                    v = omega(i)./real(solution);
                elseif fractional_step == true
                    v_att = omega(i)./real(solution_att);
                    v = omega(i)./real(solution_v);
                end

                if number_of_fluids == 0 || phi_in == 0
                    P1    = find(v==max(v));
                    v1(i) = v(P1);
                    v(P1) = 0;
                    S     = find(v==max(v));
                    v2(i) = v(S);

                    if attenuation == true
                        if fractional_step == false
                            Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                            Q2(i) = 2*abs(imag(solution(S))/real(solution(S)));
                        elseif fractional_step == true
                            P1    = find(v_att==max(v_att));
                            v_att(P1) = 0;
                            S     = find(v_att==max(v_att));
                            Q1(i) = 2*abs(imag(solution_att(P1))/real(solution_att(P1)));
                            Q2(i) = 2*abs(imag(solution_att(S))/real(solution_att(S)));
                        end
                    end
                elseif number_of_fluids == 1
                    P1    = find(v==max(v));
                    v1(i) = v(P1);
                    v(P1) = 0;
                    S     = find(v==max(v));
                    v2(i) = v(S);
                    v(S)  = 0;
                    if (fractional_step == false && length(solution) <= 4) || (fractional_step == true &&length(solution_v) <= 4)
                        v3(i) = 0;
                    else
                        P2    = find(v==max(v));
                        v3(i) = v(P2);
                    end

                    if attenuation == true
                        if fractional_step == false
                            Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                            Q2(i) = 2*abs(imag(solution(S))/real(solution(S)));
                            if length(solution) <= 4
                                Q3(i) = 0;
                            else
                                Q3(i) = 2*abs(imag(solution(P2))/real(solution(P2)));
                            end
                        elseif fractional_step == true
                            P1        = find(v_att==max(v_att));
                            Q1(i) = 2*abs(imag(solution_att(P1))/real(solution_att(P1)));
                            v_att(P1) = 0;
                            S         = find(v_att==max(v_att));
                            Q2(i) = 2*abs(imag(solution_att(S))/real(solution_att(S)));
                            v_att(S)  = 0;
                            if length(solution_att) <= 4
                                Q3(i) = 0;
                            else
                                P2        = find(v_att==max(v_att));
                                Q3(i) = 2*abs(imag(solution_att(P2))/real(solution_att(P2)));
                            end
                        end
                    end
                elseif number_of_fluids == 2
                    P1    = find(v==max(v));
                    v1(i) = v(P1);
                    v(P1) = 0;
                    S     = find(v==max(v));
                    v2(i) = v(S);
                    v(S)  = 0;
                    if (fractional_step == false && length(solution) <= 4) || (fractional_step == true && length(solution_v) <= 4)
                        v3(i) = 0;
                        v4(i) = 0;
                    else
                        P2    = find(v==max(v));
                        v3(i) = v(P2);
                        if length(solution) > 6
                            v(P2) = 0;
                            P3    = find(v==max(v));
                            v4(i) = v(P3);
                        else
                            v4(i) = 0;
                        end
                    end

                    if attenuation == true
                        if fractional_step == false
                            Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                            Q2(i) = 2*abs(imag(solution(S))/real(solution(S)));
                            if length(solution) <= 4
                                Q3(i) = 0;
                                Q4(i) = 0;
                            else
                                Q3(i) = 2*abs(imag(solution(P2))/real(solution(P2)));
                                if length(solution) > 6
                                    Q4(i) = 2*abs(imag(solution(P3))/real(solution(P3)));
                                else
                                    Q4(i) = 0;
                                end
                            end
                        elseif fractional_step == true
                            P1        = find(v_att==max(v_att));
                            Q1(i) = 2*abs(imag(solution_att(P1))/real(solution_att(P1)));
                            v_att(P1) = 0;
                            S         = find(v_att==max(v_att));
                            Q2(i) = 2*abs(imag(solution_att(S))/real(solution_att(S)));
                            v_att(S)  = 0;
                            if length(solution_att) <= 4
                                Q3(i) = 0;
                                Q4(i) = 0;
                            else
                                P2        = find(v_att==max(v_att));
                                Q3(i) = 2*abs(imag(solution_att(P2))/real(solution_att(P2)));
                                v_att(P2) = 0;
                                P3        = find(v_att==max(v_att));
                                Q4(i) = 2*abs(imag(solution_att(P3))/real(solution_att(P3)));
                            end
                        end
                    end
                end
            end

            %% Calculation of the wavespeed and the inverse quality factor for the most energetic frequency
            omega_f0 = linspace(min(omega),omega_p,2);

            if number_of_fluids == 1
                v1_fc = interp1(omega,v1,omega_f0);
                v2_fc = interp1(omega,v2,omega_f0);
                v3_fc = interp1(omega,v3,omega_f0);

                Q1_fc = interp1(omega,Q1,omega_f0);
                Q2_fc = interp1(omega,Q2,omega_f0);
                Q3_fc = interp1(omega,Q3,omega_f0);
            elseif number_of_fluids == 2
                v1_fc = interp1(omega,v1,omega_f0);
                v2_fc = interp1(omega,v2,omega_f0);
                v3_fc = interp1(omega,v3,omega_f0);
                v4_fc = interp1(omega,v4,omega_f0);

                Q1_fc = interp1(omega,Q1,omega_f0);
                Q2_fc = interp1(omega,Q2,omega_f0);
                Q3_fc = interp1(omega,Q3,omega_f0);
                Q4_fc = interp1(omega,Q4,omega_f0);
            end
            
            if number_of_fluids >= 1
                wavelength1 = v1./(omega'./2*pi);
                wavelength2 = v2./(omega'./2*pi);
                wavelength3 = v3./(omega'./2*pi);
                if number_of_fluids == 2
                    wavelength4 = v4./(omega'./2*pi);
                end
            end
            
            %% PLOTTING OUTPUT
            if plot_diag && number_of_fluids >= 1
                cd (workpath)
                run('plotting_output')
            end

            %% PRINTING OUTPUT ON SCREEN
            cd (workpath)
            run('printing_output')
            fprintf('\n');

            %% CALCULATION OF FREQUENCY STABILITY
            %s1 = min(v1) / (2.5*f0);
            %s2 = min(v2) / (2.5*f0);
            s1 = min(v1./omega') * 2*pi;
            s2 = min(v2./omega') * 2*pi;

            if number_of_fluids >= 1
                %s3 = min(v3) / (2.5*f0);
                s3 = min(v3./omega') * 2*pi;
            end
            if number_of_fluids >= 2
                %s4 = min(v4) / (2.5*f0);
                s4 = min(v4./omega') * 2*pi;
            end

            if number_of_fluids == 0
                s_frequ(u) = min([s1,s2]);
            elseif number_of_fluids == 1
                s_frequ(u) = min([s1,s2,s3]);
            elseif number_of_fluids == 2
                s_frequ(u) = min([s1,s2,s3,s4]);
            end
            fprintf('Size of an element for material %g: %g m \n',u,s_frequ(u));
            fprintf('--------------------------------------------------- \n');
            fprintf('\n');
        end
        fprintf('Size of one element to show each wave for all materials: %g m \n',min(s_frequ));
        %fprintf('Frequency stability is calculated by the given most energetic frequency in "source".\n');
        fprintf('------------------------------------------------------------------------------------- \n');
        fprintf('\n');

    end
end

if logfile
diary off
end
cd (workpath)
cd ../tools
