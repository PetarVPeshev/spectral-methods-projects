function [sw_pwr] = get_sw_power(k0, krho_sw, er, h_subs, mode)
%GET_SW_POWER Summary of this function goes here
%   k0 - wave propagation constant in free space
%   er - relative permittivity of substrate / dielectric
%   krho_tm - rho-component of the TM0 mode propagation constant
%   h_subs - height of substrate / dielectric
    eta_fs = 376.730313668;        % Free space wave impedance

    k_subs = k0 * sqrt(er);        % Propagation constant in substrate
    eta_subs = eta_fs / sqrt(er);  % Substrate wave impedance

    kz_sw_fs = -1j * sqrt( - k0 ^ 2 + krho_sw .^ 2 );
    kz_sw_subs = - 1j * sqrt( - k_subs ^ 2 + krho_sw .^ 2 );

    if strcmp(mode, 'tm')
        imped_fs = eta_fs * kz_sw_fs / k0;
        imped_subs = eta_subs * kz_sw_subs / k_subs;
    elseif strcmp(mode, 'te')
        imped_fs = eta_fs * k0 ./ kz_sw_fs;
        imped_subs = eta_subs * k_subs ./ kz_sw_subs;
    else
        error('Error. Invalid mode argument.');
    end

    imped_down = 1j * imped_subs .* tan(kz_sw_subs .* h_subs);

    disper_deriv = get_dispersion_deriv(k0, krho_sw, er, h_subs, mode);

    imped_comp = abs(imped_fs .* imped_down ./ disper_deriv) .^ 2;
    subs_sin_comp = 1 ./ ( sin(kz_sw_subs .* h_subs) .^ 2 );

    if strcmp(mode, 'tm')
        integr_subs_comp ...
            = h_subs .* ( 1 + sinc(2 * kz_sw_subs .* h_subs / pi) ) / 2;
        integr_subs = imped_comp .* ( abs(1 ./ imped_subs) .^ 2 ) ...
            .* subs_sin_comp .* integr_subs_comp / er;

        integr_fs_comp = 1 ./ ( 2 * sqrt( krho_sw .^ 2 - k0 .^ 2 ) );
        integr_fs ...
            = imped_comp .* ( abs(1 ./ imped_fs) .^ 2 ) .* integr_fs_comp;

        integr_z = (integr_subs + integr_fs) * eta_fs / k0;
    elseif strcmp(mode, 'te')
        integr_subs_comp ...
            = h_subs .* ( 1 - sinc(2 * kz_sw_subs .* h_subs / pi) ) / 2;
        integr_subs = imped_comp .* subs_sin_comp .* integr_subs_comp;

        integr_fs_comp = 1 ./ ( 2 * sqrt( krho_sw .^ 2 - k0 .^ 2 ) );
        integr_fs = imped_comp .* integr_fs_comp;

        integr_z = (integr_subs + integr_fs) / (k0 * eta_fs);
    end

    integr_phi = pi;  % Assuming elementary current source

    sw_pwr = integr_z .* integr_phi .* (krho_sw .^ 2) / (4 * pi);

    sw_pwr(isinf(sw_pwr)) = 0;
    sw_pwr(isnan(sw_pwr)) = 0;

%     %% TM0 POWER
%     kz_tm_fs = -1j * sqrt( - k0 ^ 2 + krho_tm .^ 2 );
%     kz_tm_subs = - 1j * sqrt( - k_subs ^ 2 + krho_tm .^ 2 );
% 
%     imped_tm_fs = eta_fs * kz_tm_fs / k0;  % Upward impedance
%     imped_tm_subs = eta_subs * kz_tm_subs / k_subs;
% 
%     imped_tm_down = 1j * imped_tm_subs .* tan(kz_tm_subs .* h_subs);
% 
%     disper_tm_deriv = get_dispersion_deriv(k0, krho_tm, er, h_subs, 'tm');
% 
%     integr_tm_imped_comp = abs(imped_tm_fs .* imped_tm_down ...
%         ./ disper_tm_deriv) .^ 2;
% 
%     integr_tm_subs_sin = 1 ./ ( sin(kz_tm_subs .* h_subs) .^ 2 );
%     integr_tm_subs_int = h_subs ...
%         .* ( 1 + sinc(2 * kz_tm_subs .* h_subs / pi) ) / 2;
%     integr_tm_subs = integr_tm_imped_comp ...
%         .* ( abs(1 ./ imped_tm_subs) .^ 2 ) .* integr_tm_subs_sin ...
%         .* integr_tm_subs_int / er;
% 
%     integr_tm_fs_int = 1 ./ ( 2 * sqrt( krho_tm .^ 2 - k0 .^ 2 ) );
%     integr_tm_fs = integr_tm_imped_comp ...
%         .* ( abs(1 ./ imped_tm_fs) .^ 2 ) .* integr_tm_fs_int;
% 
%     integr_z_tm = (integr_tm_subs + integr_tm_fs) * eta_fs / k0;
%     integr_phi_tm = pi;  % Assuming short dipole
% 
%     sw_tm0_pwr = integr_z_tm .* integr_phi_tm ...
%         .* (krho_tm .^ 2) / (4 * pi);
% 
%     sw_tm0_pwr(isnan(sw_tm0_pwr)) = 0;
% 
%     %% TE1 POWER
%     kz_te_fs = -1j * sqrt( - k0 ^ 2 + krho_te .^ 2 );
%     kz_te_subs = - 1j * sqrt( - k_subs ^ 2 + krho_te .^ 2 );
% 
%     imped_te_fs = eta_fs * k0 ./ kz_te_fs;  % Upward impedance
%     imped_te_subs = eta_subs * k_subs ./ kz_te_subs;
% 
%     imped_te_down = 1j * imped_te_subs .* tan(kz_te_subs .* h_subs);
% 
%     disper_te_deriv = get_dispersion_deriv(k0, krho_te, er, h_subs, 'te');
% 
%     integr_te_imped_comp = abs(imped_te_fs .* imped_te_down ...
%         ./ disper_te_deriv) .^ 2;
% 
%     integr_te_subs_sin = 1 ./ ( sin(kz_te_subs .* h_subs) .^ 2 );
%     integr_te_subs_int = h_subs ...
%         .* ( 1 - sinc(2 * kz_te_subs .* h_subs / pi) ) / 2;
%     integr_te_subs = integr_te_imped_comp .* integr_te_subs_sin ...
%         .* integr_te_subs_int;
% 
%     integr_te_fs_int = 1 ./ ( 2 * sqrt( krho_te .^ 2 - k0 .^ 2 ) );
%     integr_te_fs = integr_te_imped_comp .* integr_te_fs_int;
% 
%     integr_z_te = (integr_te_subs + integr_te_fs) / (k0 * eta_fs);
%     integr_phi_te = pi;  % Assuming short dipole
% 
%     sw_te1_pwr = integr_z_te .* integr_phi_te ...
%         .* (krho_te .^ 2) / (4 * pi);
% 
%     sw_te1_pwr(isinf(sw_te1_pwr)) = 0;
%     sw_te1_pwr(isnan(sw_te1_pwr)) = 0;
end

