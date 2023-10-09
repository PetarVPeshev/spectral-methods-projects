function eqn_result = get_dispersion_eqn(k0, krho_sw, er, h_subs, mode)
%GET_DISPERSION_EQN Summary of this function goes here
%   k0 - wave propagation constant in free space
%   er - relative permittivity of substrate / dielectric
%   krho_sw - rho-component of the surface wave mode propagation constant
%   h_subs - height of substrate / dielectric
%   mode - TM / TE surface wave mode ('te' / 'tm')
    eta_fs = 376.730313668;        % Free space wave impedance

    k_subs = k0 * sqrt(er);        % Propagation constant in substrate
    eta_subs = eta_fs / sqrt(er);  % Substrate wave impedance

    kz_sw_fs = - 1j * sqrt( - k0 ^ 2 + krho_sw .^ 2 );
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
    
    eqn_result = imped_fs + imped_down;
end
