function [v, i] = get_v_and_i_stratif(k0, krho, er, h_subs, z, ...
    mode, region)
%GET_I_AND_V_STRATIF Summary of this function goes here
%   k0 - wave propagation constant in free space
%   krho - rho-component of the propagation vector
%   er - relative permittivity of substrate / dielectric
%   h_subs - height of substrate / dielectric
%   z - elevation at which the field is evaluated
%   mode - TM / TE surface wave mode ('te' / 'tm')
%   region - medium in which the voltage and current are evaluated ('FS' / 'SUBS')
    eta_fs = 376.730313668;        % Free space wave impedance

    k_subs = k0 * sqrt(er);        % Propagation constant in substrate
    eta_subs = eta_fs / sqrt(er);  % Substrate wave impedance

    kz_fs = - 1j * sqrt( - k0 ^ 2 + krho .^ 2 );
    kz_subs = - 1j * sqrt( - k_subs ^ 2 + krho .^ 2 );
    
    if strcmp(mode, 'tm')
        imped_fs = eta_fs * kz_fs / k0;
        imped_subs = eta_subs * kz_subs / k_subs;
%         admit_subs = 1 ./ imped_subs;
    elseif strcmp(mode, 'te')
        imped_fs = eta_fs * k0 ./ kz_fs;
        imped_subs = eta_subs * k_subs ./ kz_subs;
%         admit_subs = 1 ./ imped_subs;
    else
        error('Error. Invalid mode argument.');
    end

    refl_coef = (imped_fs - imped_subs) ./ (imped_fs + imped_subs);
%     admit_in = admit_subs ...
%         .* ( imped_subs + 1j * imped_fs .* tan(kz_subs * h_subs) ) ...
%         ./ ( imped_fs + 1j * imped_subs .* tan(kz_subs * h_subs) );

    if strcmp(region, 'FS')
%         i = (imped_subs .* admit_in ./ imped_fs) ...
%             .* exp(1j * h_subs * (kz_fs - kz_subs)) ...
%             .* exp(- 1j * kz_fs .* z) .* (1 - refl_coef) ...
%             ./ ( refl_coef .* exp(- 2j * kz_subs * h_subs) + 1 );
%         v = imped_fs .* i;

        v = exp(- 1j * kz_fs .* (z - h_subs)) .* (1 + refl_coef) ...
            ./ ( refl_coef .* exp(- 1j * kz_subs * h_subs) ...
            + exp(1j * kz_subs * h_subs) );
        i = v ./ imped_fs;
    elseif strcmp(region, 'SUBS')
%         i = admit_in .* exp(- 1j * kz_subs .* z) ...
%             .* (1 + refl_coef .* exp(2j * kz_subs .* (z - h_subs))) ...
%             ./ ( refl_coef .* exp(- 2j * kz_subs * h_subs) + 1 );
%         v = imped_subs .* admit_in .* exp(- 1j * kz_subs .* z) ...
%             .* (1 - refl_coef .* exp(2j * kz_subs .* (z - h_subs))) ...
%             ./ ( refl_coef .* exp(- 2j * kz_subs * h_subs) + 1 );

        v = exp(- 1j * kz_subs .* z) ...
            .* (1 + refl_coef .* exp(2j * kz_subs .* (z - h_subs))) ...
            ./ (refl_coef .* exp(- 2j * kz_subs * h_subs) + 1);
        i = exp(- 1j * kz_subs .* z) ...
            .* (1 - refl_coef .* exp(2j * kz_subs .* (z - h_subs))) ...
            ./ ( imped_subs ...
            .* (refl_coef .* exp(- 2j * kz_subs * h_subs) + 1) );
    else
        error('Error. Invalid region argument.');
    end
end

