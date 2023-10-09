close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\figures\fields'], 'dir')
    mkdir('figures\fields');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% DIELECTRIC MATERIALS
dielectric_list;

%% WAVE & SUBSTRATE PARAMETERS
h_subs_norm = 0.25;
er_subs = 2.1;
name_subs = 'PTFE';
freq = 2.4 * 1e9;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;
wlen_subs = wlen / sqrt(er_subs);
h_subs = h_subs_norm * wlen_subs;
    
%% CYLINDRICAL COORDINATES
rho = linspace(0.05, 5, 101) * wlen_subs;
phi = linspace(0, 2 * pi, 101);
z = linspace(0, 0.5, 51) * wlen_subs;
[RHO, PHI] = meshgrid(rho, phi);

%% SURFACE WAVE MODE PROPAGATION VECTOR
krho = linspace(1, sqrt(er_subs), 1001) * k0;
krho_tm = get_krho_sw(k0, krho, er_subs, h_subs, 'tm');

%% CURRENT SOURCE FOURIER TRANSFORM
j_ft = ones( [size(PHI, 1, 2)] );
[kx_sw, ~] = get_ksw_comp(krho_tm, PHI);
dx = pi / krho_tm;
j_ft_ts = j_ft .* 2 .* cos(kx_sw .* dx / 2);

%% SURFACE WAVE FIELDS
% Two sources
[e_sw_min, h_sw_min] = get_sw_fields(k0, krho_tm, j_ft_ts, er_subs, ...
    h_subs, RHO, PHI, z, 'tm');
% Single source
[e_sw, h_sw] = get_sw_fields(k0, krho_tm, j_ft, er_subs, ...
    h_subs, RHO, PHI, z, 'tm');

%% SURFACE WAVE TOTAL FIELDS
% evaluated at rho = 2.525 * wlen (51), z = 0.12 * wlen (13)
% Two sources
et_sw_min = sqrt(abs(e_sw_min(:, 51, 13, 1)) .^ 2 ...
    + abs(e_sw_min(:, 51, 13, 2)) .^ 2 + abs(e_sw_min(:, 51, 13, 3)) .^ 2);
ht_sw_min = sqrt(abs(h_sw_min(:, 51, 13, 1)) .^ 2 ...
    + abs(h_sw_min(:, 51, 13, 2)) .^ 2 + abs(h_sw_min(:, 51, 13, 3)) .^ 2);
% Single source
et_sw = sqrt(abs(e_sw(:, 51, 13, 1)) .^ 2 ...
    + abs(e_sw(:, 51, 13, 2)) .^ 2 + abs(e_sw(:, 51, 13, 3)) .^ 2);
ht_sw = sqrt(abs(h_sw(:, 51, 13, 1)) .^ 2 ...
    + abs(h_sw(:, 51, 13, 2)) .^ 2 + abs(h_sw(:, 51, 13, 3)) .^ 2);

%% PLOT TM SURFACE WAVE E & H-FIELDS
color_styles = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", ...
"#4DBEEE"];

figure('Position', [250 250 750 400]);
plot(phi * 180 / pi, et_sw_min, 'LineWidth', 2.0, ...
    'DisplayName', 'two dipoles, d = \lambda_{sw} / 2');
hold on;
plot(phi * 180 / pi, et_sw, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'single dipole');
grid on;
xticks(0 : 30 : 360);
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta');
ylabel('|E_{SW}^{TM}|');
title(['TM Surface Wave @ j_{x} = 1, \rho = 2.52\lambda_{d}, ' ...
    'z = 0.12\lambda_{d}']);
saveas(gcf, 'figures\fields\et_sw.fig');

figure();
polarplot(phi, et_sw_min, 'LineWidth', 2.0, ...
    'DisplayName', 'two dipoles, d = \lambda_{sw} / 2');
hold on;
polarplot(phi, et_sw, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'single dipole');
legend show;
legend('location', 'bestoutside');
title(['|E_{SW}^{TM}| @ j_{x} = 1, \rho = 2.52\lambda_{d}, ' ...
    'z = 0.12\lambda_{d}']);

figure('Position', [250 250 750 400]);
plot(phi * 180 / pi, ht_sw_min, 'LineWidth', 2.0, ...
    'DisplayName', 'two dipoles, d = \lambda_{sw} / 2');
hold on;
plot(phi * 180 / pi, ht_sw, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'single dipole');
grid on;
xticks(0 : 30 : 360);
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta');
ylabel('|H_{SW}^{TM}|');
title(['TM Surface Wave @ j_{x} = 1, \rho = 2.52\lambda_{d}, ' ...
    'z = 0.12\lambda_{d}']);
saveas(gcf, 'figures\fields\ht_sw.fig');
