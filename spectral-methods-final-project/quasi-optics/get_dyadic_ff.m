function field = get_dyadic_ff(k, r, sph_grid, kz, sgf, current_ft)
%GET_DYADIC_FF This function evaluates the electric or magnetic far-field in 
%spherical coordinates for a given Fourier Transform of electric or 
%magnetic current density
%   Detailed explanation goes here
    sgf_current = permute( pagemtimes( permute(sgf, [3, 4, 1, 2]), ...
        permute(current_ft, [3, 4, 1, 2]) ), [3, 4, 1, 2]);
    field = 1j * cart2sph_vector(kz .* sgf_current, sph_grid) .* ...
        exp(-1j * k * r) ./ (2 * pi * r);
end

