close all;
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

%% SURFACE WAVE MODE & TWO SOURCE DISTANCE
krho_tm = NaN(1, length(h_subs));
krho_te = NaN(1, length(h_subs));
krho = linspace(1, sqrt(er_subs), 1001) * k0;
for h_idx = 1 : 1 : length(h_subs)
    krho_tm(h_idx) = get_krho_sw(k0, krho, er_subs, h_subs(h_idx), 'tm');
    krho_te(h_idx) = get_krho_sw(k0, krho, er_subs, h_subs(h_idx), 'te');
end

%% SURFACE WAVE POWER
sw_pwr_tm = get_sw_power(k0, krho_tm, er_subs, h_subs, 'tm');
sw_pwr_te = get_sw_power(k0, krho_te, er_subs, h_subs, 'te');

%% WAVE VECTOR
[kx, ky, kz] = get_wave_vec(1, k0, THETA, PHI);
krho = sqrt(kx .^ 2 + ky .^ 2);

rad_pwr = NaN(1, length(h_subs));
for h_idx = 1 : 1 : length(h_subs)
    %% STRATIFIED MEDIA VOLTAGE & CURRENT
    [v_te, i_te, v_tm, i_tm] = stratified_media(k0, krho, z_observ, ...
        'GroundSlab', h_subs(h_idx), er_subs);
    
    %% SPECTRAL GREEN'S FUNCTION
    sgf = spectral_gf(1, k0, kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'J');
    
    %% ELEMENTARY ELECTRIC CURRENT
    j_ft = zeros( [size(sph_grid, 1, 2), 2] );
    j_ft(:, :, 1) = 1;
    
    %% E FAR-FIELD
    e_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, j_ft);

    %% RADIATED POWER
    [~, ~, rad_pwr(h_idx)] = directivity(1, e_ff, sph_grid, r_observ);
end

%% SURFACE WAVE EFFICIENCY
sw_eff = rad_pwr ./ (rad_pwr + sw_pwr_tm + sw_pwr_te);

%% PLOT SURFACE WAVE AND RADIATED POWER
figure('Position', [250 250 750 400]);
plot(h_subs_norm, sw_pwr_tm, 'LineWidth', 2.0, 'DisplayName', 'TM0');
hold on;
plot(h_subs_norm, sw_pwr_te, 'LineWidth', 2.0, 'DisplayName', 'TE1');
hold on;
plot(h_subs_norm, rad_pwr, 'LineWidth', 2.0, 'DisplayName', 'radiated');
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P');
title(['Radiated & Surface Wave Power @ j_{x} = 1, ' name_subs ...
    ' \epsilon_{r} = ' num2str(er_subs)]);
saveas(gcf, 'figures\characterization\subs_antenna.fig');

figure('Position', [250 250 750 400]);
plot(h_subs_norm, sw_eff, 'LineWidth', 2.0, ...
    'DisplayName', '\eta_{SW}');
grid on;
xlim([0.05 max(h_subs_norm)]);
ylim([0 1]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('\eta_{SW}');
title(['Surface Wave Efficiency @ j_{x} = 1, ' name_subs ...
    ' \epsilon_{r} = ' num2str(er_subs)]);
saveas(gcf, 'figures\characterization\sw_eff.fig');
