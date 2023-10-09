function [phi_integr] = get_phi_integr(j_ft, phi, mode, source)
%GET_PHI_INTEGR Summary of this function goes here
%   j_ft - FT of current distribution
%   phi - azimuth angle
%   mode - TM / TE surface wave mode ('te' / 'tm')
%   source - current density source ('ELECT' / 'MAGN')
    dphi = phi(2) - phi(1);
    
    if strcmp(source, 'ELECT')
        if strcmp(mode, 'tm')
            phi_integr = sum( (abs(j_ft) .^ 2) .* (cos(phi) .^ 2) ) * dphi;
        elseif strcmp(mode, 'te')
            phi_integr = sum( (abs(j_ft) .^ 2) .* (sin(phi) .^ 2) ) * dphi;
        else
            error('Error. Invalid mode argument.');
        end
    elseif strcmp(source, 'MAGN')
        if strcmp(mode, 'tm')
            phi_integr = sum( (abs(j_ft) .^ 2) .* (sin(phi) .^ 2) ) * dphi;
        elseif strcmp(mode, 'te')
            phi_integr = sum( (abs(j_ft) .^ 2) .* (cos(phi) .^ 2) ) * dphi;
        else
            error('Error. Invalid mode argument.');
        end
    else
        error('Error. Invalid source argument.');
    end
end

