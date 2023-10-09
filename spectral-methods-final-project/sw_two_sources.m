close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\figures\sw_cancellation'], 'dir')
    mkdir('figures\sw_cancellation');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
h_subs_norm = linspace(0.0001, 0.5, 201);
er_subs = 2.1;
name_subs = 'PTFE';
freq = 2.4 * 1e9;
r_observ = 1;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;
wlen_subs = wlen / sqrt(er_subs);
h_subs = h_subs_norm * wlen_subs;

%% SPHERICAL COORDINATE PARAMETERS
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 501);
phi = linspace(0, 2 * pi, 2001);
[THETA, PHI] = meshgrid(theta, phi);
sph_grid = meshgrid_comb(theta, phi);

%% OBSERVATION ELEVATION
z_observ = r_observ * cos(THETA);

%% SURFACE WAVE MODES
krho = linspace(1, sqrt(er_subs), 1001) * k0;
for h_idx = length(h_subs) : -1 : 1
    krho_tm(h_idx) = get_krho_sw(k0, krho, er_subs, h_subs(h_idx), 'tm');
    krho_te(h_idx) = get_krho_sw(k0, krho, er_subs, h_subs(h_idx), 'te');
end

%% SURFACE WAVE POWER (ELEMENTARY CURRENT SOURCE)
pwr_sw_tm = get_sw_power(k0, krho_tm, er_subs, h_subs, 'tm');
pwr_sw_te = get_sw_power(k0, krho_te, er_subs, h_subs, 'te');

%% DISTANCE BETWEEN TWO SOURCES
d = pi ./ krho_tm;

%% PHI INTEGRAND OF SURFACE WAVE POWER
integr_phi_tm = NaN(1, length(h_subs));
integr_phi_te = NaN(1, length(h_subs));
for h_idx = 1 : 1 : length(h_subs)
    [kx_sw_tm, ~] = get_ksw_comp(krho_tm(h_idx), phi);
    [kx_sw_te, ~] = get_ksw_comp(krho_te(h_idx), phi);
    
    j_ft_tm_min = 2 * cos(kx_sw_tm * d(h_idx) / 2) .* ones(size(phi));
    j_ft_te_min = 2 * cos(kx_sw_te * d(h_idx) / 2) .* ones(size(phi));

    integr_phi_tm(1, h_idx) = get_phi_integr(j_ft_tm_min, phi, ...
        'tm', 'ELECT');
    integr_phi_te(1, h_idx) = get_phi_integr(j_ft_te_min, phi, ...
        'te', 'ELECT');
end

%% POWER TWO SOURCES
pwr_sw_tm_min = pwr_sw_tm .* integr_phi_tm / pi;
pwr_sw_te_min = pwr_sw_te .* integr_phi_te / pi;

%% WAVE VECTOR
[kx, ky, kz] = get_wave_vec(1, k0, THETA, PHI);
krho = sqrt(kx .^ 2 + ky .^ 2);

rad_pwr = NaN(1, length(h_subs));
rad_pwr_min = NaN(1, length(h_subs));
for h_idx = 1 : 1 : length(h_subs)
    %% STRATIFIED MEDIA VOLTAGE & CURRENT
    [v_te, i_te, v_tm, i_tm] = stratified_media(k0, krho, z_observ, ...
        'GroundSlab', h_subs(h_idx), er_subs);
    
    %% SPECTRAL GREEN'S FUNCTION
    sgf = spectral_gf(1, k0, kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'J');
    
    %% ELEMENTARY ELECTRIC CURRENT
    j_ft = zeros( [size(sph_grid, 1, 2), 2] );
    j_ft(:, :, 1) = 1;
    j_ft_min = j_ft .* 2 .* cos(kx * d(h_idx) / 2);
    
    %% E FAR-FIELD
    e_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, j_ft);
    e_ff_min = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, ...
        j_ft_min);

    %% RADIATED POWER
    [~, ~, rad_pwr(h_idx)] = directivity(1, e_ff, sph_grid, r_observ);
    [~, ~, rad_pwr_min(h_idx)] = directivity(1, e_ff_min, sph_grid, ...
        r_observ);
end

%% PLOT SURFACE WAVE POWER
figure('Position', [250 250 825 400]);
plot(h_subs_norm, pwr_sw_tm, 'Color', '#0072BD', ...
    'LineWidth', 2.0, 'DisplayName', 'TM, single source');
hold on;
plot(h_subs_norm, pwr_sw_tm_min, '--', 'Color', '#0072BD', ...
    'LineWidth', 2.0, 'DisplayName', 'TM, two sources');
hold on;
plot(h_subs_norm, pwr_sw_te, 'Color', '#D95319', ...
    'LineWidth', 2.0, 'DisplayName', 'TE, single source');
hold on;
plot(h_subs_norm, pwr_sw_te_min, '--', 'Color', '#D95319', ...
    'LineWidth', 2.0, 'DisplayName', 'TE, two sources');
hold on;
plot(h_subs_norm, rad_pwr, 'Color', '#EDB120', ...
    'LineWidth', 2.0, 'DisplayName', 'radiated, single sources');
hold on;
plot(h_subs_norm, rad_pwr_min, '--', 'Color', '#EDB120', ...
    'LineWidth', 2.0, 'DisplayName', 'radiated, two sources');
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P_{SW}');
title(['Surface Wave Power @ j_{x} = 1, \epsilon_{r,subs} = ' ...
    num2str(er_subs) ', d = \lambda_{SW} / 2']);
saveas(gcf, 'figures\sw_cancellation\dtm_x.fig');
