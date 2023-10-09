% close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\figures\characterization'], 'dir')
    mkdir('figures\characterization');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% DIELECTRIC MATERIALS
dielectric_list;

%% PARAMETERS
h_subs_norm = linspace(0.0001, 0.5, 201);
freq = 2.4 * 1e9;
r_observ = 1;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;

%% SPHERICAL COORDINATE PARAMETERS
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 501);
phi = linspace(0, 2 * pi, 2001);
[THETA, PHI] = meshgrid(theta, phi);
sph_grid = meshgrid_comb(theta, phi);

%% OBSERVATION ELEVATION
z_observ = r_observ * cos(THETA);

%% WAVE VECTOR
[kx, ky, kz] = get_wave_vec(1, k0, THETA, PHI);
krho = sqrt(kx .^ 2 + ky .^ 2);

rad_pwr = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    wlen_subs = wlen / sqrt(dielectric(subs_idx).er);
    h_subs = h_subs_norm * wlen_subs;
    
    for h_idx = 1 : 1 : length(h_subs)
        %% STRATIFIED MEDIA VOLTAGE & CURRENT
        [v_tm, i_tm] = get_v_and_i_stratif(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), z_observ, 'tm', 'FS');
        [v_te, i_te] = get_v_and_i_stratif(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), z_observ, 'te', 'FS');
        
        %% SPECTRAL GREEN'S FUNCTION
        sgf = spectral_gf(1, k0, kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'M');
        
        %% ELEMENTARY ELECTRIC CURRENT
        m_ft = zeros( [size(sph_grid, 1, 2), 2] );
        m_ft(:, :, 1) = 1;
        
        %% E FAR-FIELD
        e_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, m_ft);
    
        %% RADIATED POWER
        [~, ~, rad_pwr(subs_idx, h_idx)] = directivity(1, e_ff, ...
            sph_grid, r_observ);
    end
end

%% PLOT SURFACE WAVE AND RADIATED POWER
color_styles = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", ...
    "#4DBEEE"];

figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
% for subs_idx = 1 : 1 : 1
    plot(h_subs_norm, rad_pwr(subs_idx, :), ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', [dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P_{rad}');
title('Slot Radiated Power @ m_{x} = 1');
saveas(gcf, 'figures\characterization\slot_rad_pwr.fig');
