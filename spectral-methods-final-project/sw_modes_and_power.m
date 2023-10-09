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

%% DIELECTRIC MATERIALS
dielectric_list;

%% PARAMETERS
h_subs_norm = linspace(0.0001, 0.5, 1001);
freq = 2.4 * 1e9;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;

%% SOURCE PARAMETERS
len_dipole = wlen / 2;
wid_dipole = wlen / 20;

%% AZIMUTH ANGLE
phi = linspace(0, 2 * pi, 1001);

%% SURFACE WAVE MODES
krho_te = NaN(length(dielectric), length(h_subs_norm));
krho_tm = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    krho = linspace(1, sqrt(dielectric(subs_idx).er), 1001) * k0;
    h_subs = h_subs_norm * wlen / sqrt(dielectric(subs_idx).er);
    
    for h_idx = length(h_subs) : -1 : 1
        krho_tm(subs_idx, h_idx) = get_krho_sw(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), 'tm');
        krho_te(subs_idx, h_idx) = get_krho_sw(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), 'te');
    end
end

%% SURFACE WAVE POWER (ELEMENTARY CURRENT SOURCE)
pwr_sw_tm = NaN(length(dielectric), length(h_subs_norm));
pwr_sw_te = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    h_subs = h_subs_norm * wlen / sqrt(dielectric(subs_idx).er);
    
    pwr_sw_tm(subs_idx, :) = get_sw_power(k0, krho_tm(subs_idx, :), ...
        dielectric(subs_idx).er, h_subs, 'tm');
    pwr_sw_te(subs_idx, :) = get_sw_power(k0, krho_te(subs_idx, :), ...
        dielectric(subs_idx).er, h_subs, 'te');
end

%% SURFACE WAVE POWER (DIPOLE, NORMAL)
integr_phi_tm = NaN(length(dielectric), length(h_subs_norm));
integr_phi_te = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    for h_idx = 1 : 1 : length(h_subs)
        [kx_tm, ky_tm] = get_ksw_comp(krho_tm(subs_idx, h_idx), phi);
        [kx_te, ky_te] = get_ksw_comp(krho_te(subs_idx, h_idx), phi);

        j_ft_tm = sinc(kx_tm * len_dipole / (2 * pi)) ...
            .* sinc(ky_tm * wid_dipole / (2 * pi));
        j_ft_te = sinc(kx_te * len_dipole / (2 * pi)) ...
            .* sinc(ky_te * wid_dipole / (2 * pi));

        integr_phi_tm(subs_idx, h_idx) ...
            = get_phi_integr(j_ft_tm, phi, 'tm', 'ELECT');
        integr_phi_te(subs_idx, h_idx) ...
            = get_phi_integr(j_ft_te, phi, 'te', 'ELECT');
    end
end
pwr_sw_tm_dipole = pwr_sw_tm .* integr_phi_tm / pi;
pwr_sw_te_dipole = pwr_sw_te .* integr_phi_te / pi;

%% PLOT PROPAGATION CONSTANT
color_styles = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", ...
    "#4DBEEE"];

figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(h_subs_norm, krho_tm(subs_idx, :) / k0, ...
        '-.', 'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TM0, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
    plot(h_subs_norm, krho_te(subs_idx, :) / k0, ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TE1, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('k_{\rho} / k_{0}');
title('Normalized k_{\rho}');
saveas(gcf, 'figures\characterization\sw_modes.fig');

%% PLOT SURFACE WAVE POWER (ELEMENTARY CURRENT SOURCE)
figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(h_subs_norm, pwr_sw_tm(subs_idx, :), ...
        '-.', 'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TM0, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
    plot(h_subs_norm, pwr_sw_te(subs_idx, :), ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TE1, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P_{SW}');
title('Surface Wave Power @ j_{x} = 1');
saveas(gcf, 'figures\characterization\sw_pwr_jx.fig');

%% PLOT SURFACE WAVE POWER (DIPOLE)
figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(h_subs_norm, pwr_sw_tm_dipole(subs_idx, :), ...
        '-.', 'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TM0, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
    plot(h_subs_norm, pwr_sw_te_dipole(subs_idx, :), ...
        'Color', [color_styles(subs_idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['TE1, ' dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0.05 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P_{SW}');
title(['Surface Wave Power @ Dipole, W = ' ...
    num2str(wid_dipole / wlen) '\lambda_{0}, L = ' ...
    num2str(len_dipole / wlen) '\lambda_{0}']);
saveas(gcf, 'figures\characterization\sw_pwr_dipole.fig');
