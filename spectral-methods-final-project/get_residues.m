function [v, i] = get_residues(k0, krho_sw, er, h_subs, z, mode)
%GET_RESIDUAL Summary of this function goes here
%   k0 - wave propagation constant in free space
%   krho_sw - rho-component of the surface wave propagation vector (double)
%   er - relative permittivity of substrate / dielectric (double)
%   h_subs - height of substrate / dielectric (double)
%   z - elevation at which voltage and current is evaluated (double array)
%   mode - TM / TE surface wave mode ('te' / 'tm')
    eta_fs = 376.730313668;        % Free space wave impedance

    k_subs = k0 * sqrt(er);        % Propagation constant in substrate
    eta_subs = eta_fs / sqrt(er);  % Substrate wave impedance

    kz_sw_fs = -1j * sqrt( - k0 ^ 2 + krho_sw ^ 2 );
    kz_sw_subs = - 1j * sqrt( - k_subs ^ 2 + krho_sw ^ 2 );
    
    if strcmp(mode, 'tm')
        imped_fs = eta_fs * kz_sw_fs / k0;
        imped_subs = eta_subs * kz_sw_subs / k_subs;
    elseif strcmp(mode, 'te')
        imped_fs = eta_fs * k0 / kz_sw_fs;
        imped_subs = eta_subs * k_subs / kz_sw_subs;
    else
        error('Error. Invalid mode argument.');
    end
    
    imped_down = 1j * imped_subs * tan(kz_sw_subs * h_subs);
    disper_deriv = get_dispersion_deriv(k0, krho_sw, er, h_subs, mode);
    imped = imped_fs * imped_down / disper_deriv;

    v = NaN( [size(z)] );
    i = NaN( [size(z)] );

    v(z <= h_subs) = imped * sin(kz_sw_subs * z(z <= h_subs)) ...
        / sin(kz_sw_subs * h_subs);
    i(z <= h_subs) = (imped / imped_subs) * 1j ...
        * cos(kz_sw_subs * z(z <= h_subs)) / sin(kz_sw_subs * h_subs);

    v(z > h_subs) = imped * exp(1j * kz_sw_fs * h_subs) ...
        * exp(- 1j * kz_sw_fs * z(z > h_subs));
    i(z > h_subs) = v(z > h_subs) / imped_fs;
end

