function SGF = get_dyadic_sgf(er, k, kx, ky, kz, field, current)
%GET_DYADIC_SGF This function calculates the Spectrum-Green-Functions dyad
%   Detailed explanation goes here
    eta_medium = 376.730313668 / sqrt(er);

    const = - eta_medium ./ ( 2 * k * kz );
    
    SGF = NaN( [size(kx, 1, 2), 3, 3] );
    if strcmp(field, 'E')
        if strcmp(current, 'J')
            k_comp = NaN( [size(kx, 1, 2), 3] );
            k_comp(:, :, 1) = kx;
            k_comp(:, :, 2) = ky;
            k_comp(:, :, 3) = kz;

            for row = 1 : 1 : 3
                for col = 1 : 1 : 3
                    if row == col
                        SGF(:, :, row, col) = const .* ...
                            ( k^2 - k_comp(:, :, row).^2 );
                    else
                        % only for z > 0
                        SGF(:, :, row, col) = const .* ...
                            ( - k_comp(:, :, row) .* ...
                            k_comp(:, :, col) );
                    end
                end
            end

        elseif strcmp(current, 'M')
            const = - 1 ./ ( 2 * kz );

            SGF(:, :, 1, 1) = 0;
            SGF(:, :, 1, 2) = 1j * const .* kz;
            SGF(:, :, 1, 3) = - 1j * const .* ky;

            SGF(:, :, 2, 1) = - 1j * const .* kz;
            SGF(:, :, 2, 2) = 0;
            SGF(:, :, 2, 3) = 1j * const .* kx;

            SGF(:, :, 3, 1) = 1j * const .* ky;
            SGF(:, :, 3, 2) = - 1j * const .* kx;
            SGF(:, :, 3, 3) = 0;
            
        else
            error('Error. Invalid argument.');
        end
    else
        error('Not implemented');
    end
end