function [ksw_x, ksw_y] = get_ksw_comp(krho, phi)
%GET_KSW_COMP Summary of this function goes here
%   ksw - rho-component of the surface wave propagation vector
%   phi - azimuth angle
    ksw_x = krho * cos(phi);
    ksw_y = krho * sin(phi);
end

