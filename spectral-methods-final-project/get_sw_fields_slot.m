function [e_sw, h_sw] = get_sw_fields_slot(k0, krho_sw, m_ft, er, ...
    h_subs, rho, phi, z, mode)
%GET_SW_FIELDS_SLOT Summary of this function goes here
%   Detailed explanation goes here
    eta_fs = 376.730313668;        % Free space wave impedance

    C = 1j * sqrt(krho_sw / (2 * pi)) * exp(1j * pi / 4);

    [v, i] = get_residues(k0, krho_sw, er, h_subs, z, mode);

    e_sw = zeros( [size(rho, 1, 2), length(z), 3] );
    h_sw = zeros( [size(rho, 1, 2), length(z), 3] );

    v = permute( repmat(v, [1 1 size(rho, 1, 2)]), [3 4 2 1] );
    i = permute( repmat(i, [1 1 size(rho, 1, 2)]), [3 4 2 1] );

    m_ft = repmat(m_ft, [1 1 length(z)]);

    if strcmp(mode, 'tm')
        cylin_wave_term = C .* sin(phi) .* exp(- 1j * krho_sw * rho) ...
            ./ sqrt(rho);
        cylin_wave_term = repmat(cylin_wave_term, [1 1 length(z)]);

        e_sw(:, :, :, 1) = v .* m_ft .* cylin_wave_term;
        e_sw(:, :, :, 3) = - (eta_fs * krho_sw / k0) .* i .* m_ft ...
            .* cylin_wave_term;
        e_sw(:, :, z <= h_subs, 3) = e_sw(:, :, z <= h_subs, 3) / er;

        h_sw(:, :, :, 2) = i .* m_ft .* cylin_wave_term;
    elseif strcmp(mode, 'te')
        cylin_wave_term = C .* cos(phi) .* exp(- 1j * krho_sw * rho) ...
            ./ sqrt(rho);
        cylin_wave_term = repmat(cylin_wave_term, [1 1 length(z)]);

        e_sw(:, :, :, 2) = - v .* m_ft .* cylin_wave_term;

        h_sw(:, :, :, 1) = i .* m_ft .* cylin_wave_term;
        h_sw(:, :, :, 3) = - (krho_sw / (eta_fs * k0)) .* i .* m_ft ...
            .* cylin_wave_term;
    else
        error('Error. Invalid argument.');
    end
end

