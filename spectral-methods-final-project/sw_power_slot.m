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

%% PLOT SURFACE WAVE POWER
color_styles = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", ...
    "#4DBEEE"];

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
xlim([0 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('P_{SW}');
title('Surface Wave Power @ m_{x} = 1');
saveas(gcf, 'figures\characterization\sw_pwr_mx.fig');
