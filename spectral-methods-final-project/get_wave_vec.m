function [kx, ky, kz] = get_wave_vec(er, k0, theta, phi)
%GET_WAVE_VEC This function calculates the wave vector components
%   er - permittivity of propagation medium
%   k0 - propagation constant in free space
%   theta - elevation angle
%   phi - azimuth angle
    k = k0 * sqrt(er);

    kx = k * sin(theta) .* cos(phi);
    ky = k * sin(theta) .* sin(phi);

    kz = - 1j * sqrt( - k .^ 2 + kx .^ 2 + ky .^ 2 );
end