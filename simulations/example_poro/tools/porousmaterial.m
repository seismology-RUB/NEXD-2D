% This script reads the parameters from the file 'matpropporo'
% version='$Rev: 21 $ ($Date: 2018-07-12 10:05:13 +0200 (Do, 12 Jul 2018) $, $Author: Janis Heuel, Marc S. Boxberg $)'


cd (workpath)
cd ../cubit

figure_number       = figure_number+u;
convert_matpropporo = str2num(cell2mat(matpropporo(u+1)));
plot_diag           = save_plot_diag;

used_model = convert_matpropporo(1); % Number of the used model
if ~ismember(used_model,[1 2 3 4 5 6])
    fprintf('Unknown model \n');
    cd (workpath)
    cd ../tools
    return
end

%% READING THE PARAMETERS FOR POROUS MATERIALS
rhos_in = convert_matpropporo(2);

if ismember(used_model,[1,2]) == true
    lambda_u_in = convert_matpropporo(3);
elseif ismember(used_model,[3,4,5,6]) == true
    K_d_in = convert_matpropporo(3);
end

my_in    = convert_matpropporo(4);
if poroelastic == true
    phi_in = convert_matpropporo(5);
elseif poroelastic == false
    phi_in = 0;
end
kappa_in = convert_matpropporo(6);
b_in     = convert_matpropporo(7);

if calculate_tortuosity == true
    r = convert_matpropporo(8);
    if r > 1
        error('The value of r is greater than 1. The script stops here.')
        cd (workpath)
    end
    inv_T_in = (1-r*(1-1/phi_in)).^(-1);
elseif calculate_tortuosity == false
    inv_T_in = convert_matpropporo(8);
end

if ismember(used_model,[1,2,5,6]) == true
    inv_N_in = convert_matpropporo(9);
elseif ismember(used_model,[3,4]) == true
    K_s_in = convert_matpropporo(9);
end

rho1_in = convert_matpropporo(10);
S1_in   = convert_matpropporo(11);
K1_in   = convert_matpropporo(12);
eta1_in = convert_matpropporo(13);
rho2_in = convert_matpropporo(14);
S2_in   = convert_matpropporo(15);

if number_of_fluids == 1
   if (S1_in == 0 && S2_in == 0) || phi_in == 0
       number_of_fluids = 0;
   end
elseif number_of_fluids == 2
   if (S1_in == 0 && S2_in == 0) || phi_in == 0
       number_of_fluids = 0;
   elseif (S1_in == 1 && S2_in == 0) || (S1_in == 0 && S2_in == 1)
       number_of_fluids = 1;
   end
end

K2_in   = convert_matpropporo(16);
eta2_in = convert_matpropporo(17);

if ismember(used_model,[1,3,5]) == true
    fitting_n_in   = convert_matpropporo(18);
    fitting_chi_in = convert_matpropporo(19);
elseif ismember(used_model,[2,4,6]) == true
    p_b_in         = convert_matpropporo(18);
    lambda_BC_in   = convert_matpropporo(19);
    if lambda_BC_in == 0
        phi_in = 0;
    end
elseif ismember(used_model,[7,8,9]) == true
    A_in           = convert_matpropporo(18);
end

%Read residial saturation  MB MB
Sr1_in   = convert_matpropporo(20);
Sr2_in   = convert_matpropporo(21);
%Test saturations
if Sr1_in > S1_in
    cd (workpath)
    error('The residual saturation of fluid 1 is larger than the actual saturation of fluid 1. Please change the entries in matpropporo.')
elseif Sr2_in > S2_in
    cd (workpath)
    error('The residual saturation of fluid 2 is larger than the actual saturation of fluid 2. Please change the entries in matpropporo.')
end

g = 9.81; % gravity acceleration [m/s^2]

if number_of_fluids == 0 | ~plot_diag
    k = 2;
elseif number_of_fluids == 1
    if eta2_in == 0
        k = 2;
    elseif eta2_in ~= 0
        k = k_default;
    end
elseif number_of_fluids == 2
    if eta1_in == 0 && eta2_in == 0
        k = 2;
    elseif eta1_in ~= 0 || eta2_in ~= 0
        k = k_default; % Creates the mesh. If k has high values the calculation needs more time
    end
end

% Generating the vector of the angular frequency
exponenten  = linspace(0, log10(omega_max),k);
omega       = 10.^(exponenten);  % Vector of the angular frequency (in radians per second)
Frequency   = omega./(2*pi);     % Frequency in Hz

%% CREATING VECTORS OF ALL PARAMETERS
rhos = linspace(rhos_in ,rhos_in ,k);                   % density solid
if ismember(used_model,[1,2]) == true
    lambda_u = linspace(lambda_u_in, lambda_u_in, k);   % first Lamé parameter (undrained)
elseif ismember(used_model,[3,4,5,6]) == true
    K_d = linspace(K_d_in, K_d_in, k);                  % bulk modulus of skeleton (drained)
end
my    = linspace(my_in, my_in, k);                      % shear modulus (second Lamé parameter)
phi   = linspace(phi_in, phi_in, k);                    % porosity
kappa = linspace(kappa_in, kappa_in, k);                % permeability
b     = linspace(b_in, b_in, k);                        % Biot coefficient
inv_T = linspace(inv_T_in, inv_T_in, k);                % Inverse Tortuosity
if ismember(used_model,[1,2,5,6]) == true
    inv_N = linspace(inv_N_in, inv_N_in, k);            % Inverse Biot modulus
elseif ismember(used_model,[3,4]) == true
    K_s = linspace(K_s_in, K_s_in, k);                  % bulk modlus of solid grain materual
end
rho1 = linspace(rho1_in, rho1_in, k);                   % density fluid 1
S1   = linspace(S1_in, S1_in, k);                       % saturation fluid 1
K1   = linspace(K1_in, K1_in, k);                       % bulk modolus fluid 1
eta1 = linspace(eta1_in, eta1_in, k);                   % viscosity fluid 1
rho2 = linspace(rho2_in, rho2_in, k);                   % density fluid 2
S2   = linspace(S2_in, S2_in, k);                       % saturation fluid 2
K2   = linspace(K2_in, K2_in, k);                       % bulk modulus fluid 2
eta2 = linspace(eta2_in, eta2_in, k);                   % viscosity fluid 2
if ismember(used_model,[1,3,5]) == true                 % fitting parameter for van Genuchten model n, chi
    fitting_n   = linspace(fitting_n_in, fitting_n_in, k);
    fitting_chi = linspace(fitting_chi_in, fitting_chi_in, k);
elseif ismember(used_model,[2,4,6]) == true
    p_b       = linspace(p_b_in, p_b_in, k);             % bubbling pressure
    lambda_BC = linspace(lambda_BC_in, lambda_BC_in, k); % fitting parameter for Brooks & Corey (1964) model
elseif ismember(used_model,[7,8,9]) == true
    % Douglas Jr. et al.
    A         = linspace(A_in, A_in, k);             % capillary pressure coefficient
end
Sr1         = linspace(Sr1_in, Sr1_in, k);           % residual saturation fluid 1
Sr2         = linspace(Sr2_in, Sr2_in, k);           % residual saturation fluid 2
