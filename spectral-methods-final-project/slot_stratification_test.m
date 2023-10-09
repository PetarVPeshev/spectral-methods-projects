close all;
clear;
clc;

addpath('../../spectral-methods-library');
addpath([pwd() '/quasi-optics']);
c = physconst('LightSpeed');

%% WAVE & SUBSTRATE PARAMETERS
h_subs_norm = 0.25;
er_subs = 1;
freq = 2.4 * 1e9;
r_observ = 1;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;
h_subs = h_subs_norm * wlen;

%% SPHERICAL COORDINATE PARAMETERS
theta = linspace(eps, pi / 2 - eps, 1001);
phi = linspace(0, 2 * pi, 4001);
[THETA, PHI] = meshgrid(theta, phi);
sph_grid = meshgrid_comb(theta, phi);

%% OBSERVATION ELEVATION
z_observ = r_observ * cos(THETA);

%% WAVE VECTOR
[kx, ky, kz] = get_wave_vec(1, k0, THETA, PHI);
krho = sqrt(kx .^ 2 + ky .^ 2);

%% STRATIFIED MEDIA VOLTAGE & CURRENT
[v_tm, i_tm, v_te, i_te] = get_stratif(k0, krho, er_subs, h_subs, z_observ);

%% SPECTRAL GREEN'S FUNCTION
sgf = spectral_gf(1, k0, kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'M');
sgf_dyad = get_dyadic_sgf(1, k0, kx, ky, kz, 'E', 'M');

%% ELEMENTARY MAGNETIC CURRENT
m_ft = zeros( [size(sph_grid, 1, 2), 3] );
m_ft(:, :, 1) = 1;

%% E FAR-FIELD
% Field with stratified media SGF
e_ff = farfield(k0, r_observ, sph_grid, kz, z_observ, ...
    sgf, m_ft(:, :, 1 : 2));
et_ff = total_field(e_ff);
% Field with quasi-optics SGF
e_dyad_ff = get_dyadic_ff(k0, r_observ, sph_grid, kz, sgf_dyad, m_ft);
et_dyad_ff = total_field(e_dyad_ff);

%% PLOT E FAR-FIELD
theta = NaN(1, 2001);
theta(1 : 1000) = - fliplr(THETA(1, 2 : 1001)) * 180 / pi;
theta(1001 : 2001) = THETA(1, :) * 180 / pi;

% Stratified media field
et_ff_0deg = NaN(1, 2001);
et_ff_0deg(1 : 1000) = fliplr(et_ff(2001, 2 : 1001));
et_ff_0deg(1001 : 2001) = et_ff(1, :);
% Quasi-optics field
et_dyad_ff_0deg = NaN(1, 2001);
et_dyad_ff_0deg(1 : 1000) = fliplr(et_dyad_ff(2001, 2 : 1001));
et_dyad_ff_0deg(1001 : 2001) = et_dyad_ff(1, :);

% Stratified media field
et_ff_45deg = NaN(1, 2001);
et_ff_45deg(1 : 1000) = fliplr(et_ff(2501, 2 : 1001));
et_ff_45deg(1001 : 2001) = et_ff(501, :);
et_ff_45deg(2001) = NaN;
% Quasi-optics field
et_dyad_ff_45deg = NaN(1, 2001);
et_dyad_ff_45deg(1 : 1000) = fliplr(et_dyad_ff(2501, 2 : 1001));
et_dyad_ff_45deg(1001 : 2001) = et_dyad_ff(501, :);
et_dyad_ff_45deg(2001) = NaN;

% Stratified media field
et_ff_90deg = NaN(1, 2001);
et_ff_90deg(1 : 1000) = fliplr(et_ff(3001, 2 : 1001));
et_ff_90deg(1001 : 2001) = et_ff(1001, :);
% Quasi-optics field
et_dyad_ff_90deg = NaN(1, 2001);
et_dyad_ff_90deg(1 : 1000) = fliplr(et_dyad_ff(3001, 2 : 1001));
et_dyad_ff_90deg(1001 : 2001) = et_dyad_ff(1001, :);

figure('Position', [250 250 750 400]);
plot(theta, et_ff_0deg ./ max(et_ff_0deg), ...
    'LineWidth', 2.0, 'DisplayName', 'stratified media SGF');
hold on;
plot(theta, et_dyad_ff_0deg ./ max(et_dyad_ff_0deg), ...
    '--', 'LineWidth', 2.0, 'DisplayName', 'quasi-optics SGF');
grid on;
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / |E_{max}|');
title('\phi = 0 deg, R = 1 m');

figure('Position', [250 250 750 400]);
plot(theta, et_ff_45deg ./ max(et_ff_45deg), ...
    'LineWidth', 2.0, 'DisplayName', 'stratified media SGF');
hold on;
plot(theta, et_dyad_ff_45deg ./ max(et_dyad_ff_45deg), ...
    '--', 'LineWidth', 2.0, 'DisplayName', 'quasi-optics SGF');
grid on;
xlim([-90 90]);
ylim([0 1]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / |E_{max}|');
title('\phi = 45 deg, R = 1 m');

figure('Position', [250 250 750 400]);
plot(theta, et_ff_90deg ./ max(et_ff_90deg), ...
    'LineWidth', 2.0, 'DisplayName', 'stratified media SGF');
hold on;
plot(theta, et_dyad_ff_90deg ./ max(et_dyad_ff_90deg), ...
    '--', 'LineWidth', 2.0, 'DisplayName', 'quasi-optics SGF');
grid on;
xlim([-90 90]);
ylim([0 1]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / |E_{max}|');
title('\phi = 90 deg, R = 1 m');
