close all;
clear;
clc;

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% DIELECTRIC MATERIALS
dielectric_list;

%% WAVE & SUBSTRATE PARAMETERS
h_subs_norm = linspace(0.0001, 0.5, 1001);
er_subs = 9.8;
freq = 2.4 * 1e9;

%% DEPENDENT PARAMETERS
wlen = c / freq;
k0 = 2 * pi / wlen;

%% SURFACE WAVE MODES FOR DIPOLE ABOVE DIELECTRIC
krho_te = NaN(length(dielectric), length(h_subs_norm));
krho_tm = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    krho = linspace(1, sqrt(dielectric(subs_idx).er), 1001) * k0;
    h_subs = h_subs_norm * wlen / sqrt(dielectric(subs_idx).er);
    
    for h_idx = length(h_subs) : -1 : 1
        [krho_te(subs_idx, h_idx), krho_tm(subs_idx, h_idx)] ...
            = find_krho(k0, krho, 'GroundSlab', h_subs(h_idx), ...
            dielectric(subs_idx).er);
    end
end

%% SURFACE WAVE MODES FOR SLOT BELOW DIELECTRIC
krho_te_slot = NaN(length(dielectric), length(h_subs_norm));
krho_tm_slot = NaN(length(dielectric), length(h_subs_norm));
for subs_idx = 1 : 1 : length(dielectric)
    krho = linspace(1, sqrt(dielectric(subs_idx).er), 1001) * k0;
    h_subs = h_subs_norm * wlen / sqrt(dielectric(subs_idx).er);
    
    for h_idx = length(h_subs) : -1 : 1
        krho_tm_slot(subs_idx, h_idx) = get_krho_sw(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), 'tm');
        krho_te_slot(subs_idx, h_idx) = get_krho_sw(k0, krho, ...
            dielectric(subs_idx).er, h_subs(h_idx), 'te');
    end
end

%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(h_subs_norm, krho_tm(subs_idx, :) / k0, ...
        'LineWidth', 2.0, 'DisplayName', ['dipole, ' ...
        dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
    plot(h_subs_norm, krho_tm_slot(subs_idx, :) / k0, ...
        '--', 'LineWidth', 2.0, 'DisplayName', ['slot, ' ...
        dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('k_{\rho} / k_{0}');
title('Normalized TM0 k_{\rho}');

figure('Position', [250 250 750 400]);
for subs_idx = 1 : 1 : length(dielectric)
    plot(h_subs_norm, krho_te(subs_idx, :) / k0, ...
        'LineWidth', 2.0, 'DisplayName', ['dipole, ' ...
        dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
    plot(h_subs_norm, krho_te_slot(subs_idx, :) / k0, ...
        '--', 'LineWidth', 2.0, 'DisplayName', ['slot, ' ...
        dielectric(subs_idx).name ' (e_{r} = ' ...
        num2str(dielectric(subs_idx).er) ')']);
    hold on;
end
grid on;
xlim([0 max(h_subs_norm)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('k_{\rho} / k_{0}');
title('Normalized TE1 k_{\rho}');
