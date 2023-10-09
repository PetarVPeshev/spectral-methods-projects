close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\figures\rad_pattern'], 'dir')
    mkdir('figures\rad_pattern');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
h_subs_norm = 0.25;
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
theta = linspace(eps, pi / 2 - eps, 1001);
phi = linspace(0, 2 * pi, 4001);
[THETA, PHI] = meshgrid(theta, phi);
sph_grid = meshgrid_comb(theta, phi);

%% OBSERVATION ELEVATION
z_observ = r_observ * cos(THETA);

%% SURFACE WAVE MODE & TWO SOURCE DISTANCE
krho = linspace(1, sqrt(er_subs), 1001) * k0;
krho_tm = get_krho_sw(k0, krho, er_subs, h_subs, 'tm');
dx = pi / krho_tm;

%% WAVE VECTOR
[kx, ky, kz] = get_wave_vec(1, k0, THETA, PHI);
krho = sqrt(kx .^ 2 + ky .^ 2);

%% STRATIFIED MEDIA VOLTAGE & CURRENT
[v_te, i_te, v_tm, i_tm] = stratified_media(k0, krho, z_observ, ...
    'GroundSlab', h_subs, er_subs);

%% SPECTRAL GREEN'S FUNCTION
sgf = spectral_gf(1, k0, kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'J');

%% ELEMENTARY ELECTRIC CURRENT
% Single source
j_ft = zeros( [size(sph_grid, 1, 2), 2] );
j_ft(:, :, 1) = 1;
% Two sources, cancelling TM0 wave
[kx_sw, ~] = get_ksw_comp(krho_tm, PHI);
j_ts_ft = 2 * cos(kx_sw * dx / 2) .* j_ft;

%% E FAR-FIELD
% Single source
e_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, j_ft);
et_ff = total_field(e_ff);
% Two sources
e_ts_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, sgf, j_ts_ft);
et_ts_ff = total_field(e_ts_ff);

%% PLOT E FAR-FIELD
theta = NaN(1, 2001);
theta(1 : 1000) = - fliplr(THETA(1, 2 : 1001)) * 180 / pi;
theta(1001 : 2001) = THETA(1, :) * 180 / pi;

e_ff_max = max(et_ff, [], 'all');
if e_ff_max < max(et_ts_ff, [], 'all')
    e_ff_max = max(et_ts_ff, [], 'all');
end

% Single source
et_ff_0deg = NaN(1, 2001);
et_ff_0deg(1 : 1000) = fliplr(et_ff(2001, 2 : 1001));
et_ff_0deg(1001 : 2001) = et_ff(1, :);
% Two sources
et_ts_ff_0deg = NaN(1, 2001);
et_ts_ff_0deg(1 : 1000) = fliplr(et_ts_ff(2001, 2 : 1001));
et_ts_ff_0deg(1001 : 2001) = et_ts_ff(1, :);

% Single source
et_ff_45deg = NaN(1, 2001);
et_ff_45deg(1 : 1000) = fliplr(et_ff(2501, 2 : 1001));
et_ff_45deg(1001 : 2001) = et_ff(501, :);
et_ff_45deg(2001) = NaN;
% Two sources
et_ts_ff_45deg = NaN(1, 2001);
et_ts_ff_45deg(1 : 1000) = fliplr(et_ts_ff(2501, 2 : 1001));
et_ts_ff_45deg(1001 : 2001) = et_ts_ff(501, :);

% Single source
et_ff_90deg = NaN(1, 2001);
et_ff_90deg(1 : 1000) = fliplr(et_ff(3001, 2 : 1001));
et_ff_90deg(1001 : 2001) = et_ff(1001, :);
% Two sources
et_ts_ff_90deg = NaN(1, 2001);
et_ts_ff_90deg(1 : 1000) = fliplr(et_ts_ff(3001, 2 : 1001));
et_ts_ff_90deg(1001 : 2001) = et_ts_ff(1001, :);

figure('Position', [250 250 750 400]);
plot(theta, 10 * log10(et_ff_0deg ./ max(e_ff_max)), ...
    'LineWidth', 2.0, 'DisplayName', 'single source');
hold on;
plot(theta, 10 * log10(et_ts_ff_0deg ./ max(e_ff_max)), ...
    '--', 'LineWidth', 2.0, 'DisplayName', ['two sources, d_{x} = ' ...
    '\lambda_{SW,TM} / 2']);
grid on;
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta [deg]');
ylabel('|E| / |E_{max}| [dB]');
title('Radiation Pattern @ j_{x} = 1, \phi = 0 deg, R = 1 m');
saveas(gcf, 'figures\rad_pattern\ej_phi_0deg.fig');

figure('Position', [250 250 750 400]);
plot(theta, 10 * log10(et_ff_45deg ./ max(e_ff_max)), ...
    'LineWidth', 2.0, 'DisplayName', 'single source');
hold on;
plot(theta, 10 * log10(et_ts_ff_45deg ./ max(e_ff_max)), ...
    '--', 'LineWidth', 2.0, 'DisplayName', ['two sources, d_{x} = ' ...
    '\lambda_{SW,TM} / 2']);
grid on;
xlim([-90 90]);
ylim([-20 0]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta [deg]');
ylabel('|E| / |E_{max}| [dB]');
title('Radiation Pattern @ j_{x} = 1, \phi = 45 deg, R = 1 m');
saveas(gcf, 'figures\rad_pattern\ej_phi_45deg.fig');

figure('Position', [250 250 750 400]);
plot(theta, 10 * log10(et_ff_90deg ./ max(e_ff_max)), ...
    'LineWidth', 2.0, 'DisplayName', 'single source');
hold on;
plot(theta, 10 * log10(et_ts_ff_90deg ./ max(e_ff_max)), ...
    '--', 'LineWidth', 2.0, 'DisplayName', ['two sources, d_{x} = ' ...
    '\lambda_{SW,TM} / 2']);
grid on;
annotation('doublearrow', [0.365, 0.365], [0.925, 0.805], 'Color', 'red');
text(5, -1.5, '- 3 dB', 'Color', 'red');
xlim([-90 90]);
ylim([-20 0]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta [deg]');
ylabel('|E| / |E_{max}| [dB]');
title('Radiation Pattern @ j_{x} = 1, \phi = 90 deg, R = 1 m');
saveas(gcf, 'figures\rad_pattern\ej_phi_90deg.fig');
