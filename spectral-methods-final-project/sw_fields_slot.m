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
freq = 2.4 * 1e9;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;

et_sw = NaN(101, length(dielectric));
ht_sw = NaN(101, length(dielectric));
for subs_idx = 1 : 1 : length(dielectric)
    wlen_subs = wlen / dielectric(subs_idx).er;
    
    %% CYLINDRICAL COORDINATES
    rho = linspace(0.05, 5, 101) * wlen_subs;
    phi = linspace(0, 2 * pi, 101);
    z = linspace(0, 0.5, 51) * wlen_subs;
    [RHO, PHI] = meshgrid(rho, phi);

    %% SURFACE WAVE MODE PROPAGATION VECTOR
    krho = linspace(1, sqrt(dielectric(subs_idx).er), 1001) * k0;
    h_subs = h_subs_norm * wlen_subs;
    krho_tm = get_krho_sw(k0, krho, dielectric(subs_idx).er, h_subs, 'tm');

    %% CURRENT SOURCE FOURIER TRANSFORM
    m_ft = ones( [size(PHI, 1, 2)] );
    [~, ky_sw] = get_ksw_comp(krho_tm, PHI);
    dy = pi / krho_tm;
    m_ft = m_ft .* 2 .* cos(ky_sw .* dy / 2);

    %% SURFACE WAVE FIELDS
    [e_sw, h_sw] = get_sw_fields_slot(k0, krho_tm, m_ft, ...
        dielectric(subs_idx).er, h_subs, RHO, PHI, z, 'tm');

    %% SURFACE WAVE TOTAL FIELDS
    % evaluated at rho = 2.525 * wlen (51), z = 0.12 * wlen (13)
    et_sw(:, subs_idx) = sqrt(abs(e_sw(:, 51, 13, 1)) .^ 2 ...
        + abs(e_sw(:, 51, 13, 2)) .^ 2 + abs(e_sw(:, 51, 13, 3)) .^ 2);
    ht_sw(:, subs_idx) = sqrt(abs(h_sw(:, 51, 13, 1)) .^ 2 ...
        + abs(h_sw(:, 51, 13, 2)) .^ 2 + abs(h_sw(:, 51, 13, 3)) .^ 2);
end

%% PLOT TM SURFACE WAVE E & H-FIELDS
color_styles = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", ...
"#4DBEEE"];

figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(phi * 180 / pi, et_sw(:, subs_idx), ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', [dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xticks(0 : 15 : 360);
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta');
ylabel('|E_{SW}^{TM}|');
title(['TM Surface Wave Total E-Field @ h = 0.25\lambda_{d}, ' ...
    '\rho = 2.525\lambda_{d}, z = 0.12\lambda_{d}']);
saveas(gcf, 'figures\fields\et_sw_mx.fig');

figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(phi * 180 / pi, ht_sw(:, subs_idx), ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', [dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xticks(0 : 15 : 360);
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta');
ylabel('|H_{SW}^{TM}|');
title(['TM Surface Wave Total H-Field @ h = 0.25\lambda_{d}, ' ...
    '\rho = 2.525\lambda_{d}, z = 0.12\lambda_{d}']);
saveas(gcf, 'figures\fields\ht_sw_mx.fig');
