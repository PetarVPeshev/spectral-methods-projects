function deriv_result = get_dispersion_deriv(k0, krho_sw, er, h_subs, mode)
%GET_DISPERSION_DERIV Summary of this function goes here
%   k0 - wave propagation constant in free space
%   krho_sw - rho-component of the surface wave mode propagation vector
%   er - relative permittivity of substrate / dielectric
%   h_subs - height of substrate / dielectric
%   mode - TM / TE surface wave mode ('te' / 'tm')
    k_delta = k0 / 500;
    
    krho_sw_high = krho_sw + k_delta / 2;
    disper_high = get_dispersion_eqn(k0, krho_sw_high, er, h_subs, mode);

    krho_sw_low = krho_sw - k_delta / 2;
    disper_low = get_dispersion_eqn(k0, krho_sw_low, er, h_subs, mode);

    deriv_result = ( disper_high - disper_low ) / k_delta;
end
