function krho_sw = get_krho_sw(k0, krho, er, h_subs, mode)
%GET_KRHO_SW Summary of this function goes here
%   Detailed explanation goes here
    disper_eqn = get_dispersion_eqn(k0, krho, er, h_subs, mode);

    % Guess krho
    [~, peak_idx] = findpeaks( abs(1 ./ disper_eqn) );
    krho_sw = krho(peak_idx);

    % Step k
    delta_k = k0 / 500;
        
    % Newton's method for SW krho
    krho_sw_prev = 0;
    if isempty(krho_sw)
        krho_sw = k0;
        krho_sw_prev = k0;
    end

    if length(krho_sw) ~= 1
        krho_sw = krho_sw(end);
        warning('Several solutions to dispersion equation are found.');
    end

    while abs(krho_sw - krho_sw_prev) > 0.00001
        disper_eqn = get_dispersion_eqn(k0, krho_sw, er, h_subs, mode);
    
        disper_eqn_1 = get_dispersion_eqn(k0, krho_sw + delta_k / 2, ...
            er, h_subs, mode);
        disper_eqn_2 = get_dispersion_eqn(k0, krho_sw - delta_k / 2, ...
            er, h_subs, mode);
    
        disper_eqn_deriv = (disper_eqn_1 - disper_eqn_2) / delta_k;
        krho_sw_prev = krho_sw;
        krho_sw = krho_sw - disper_eqn / disper_eqn_deriv;
        
        if imag(krho_sw) ~= 0
            krho_sw = k0;
            break;
        end
    end
end

