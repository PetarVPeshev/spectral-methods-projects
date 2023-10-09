close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
bottom_medium.h = 5e-3;
dielectric.er = [2.5 6 12];
wave.f  = 30 * 1e9;
N = 10001;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% WAVE PROPAGATION CONSTANT
k = wave.k0 * sqrt(dielectric.er);
    
%% Z COORDINATE
z = ones(1, N) * bottom_medium.h;

%% WAVE VECTOR COMPONENTS
kx = linspace(0, 2 * wave.k0, N);
ky = zeros( size(kx, 1, 2) );
krho = sqrt(kx .^ 2 + ky .^ 2);

SGF = NaN(1, N, 3, 2, length(dielectric.er));
for idx = 1 : 1  : length(dielectric.er)
    %% STRATIFIED MEDIA VOLTAGES AND CURRENTS
    [v_te, i_te, v_tm, i_tm] = stratified_media(wave.k0, krho, z, ...
        'SemiInfiniteSuperstrate', bottom_medium.h, dielectric.er(idx));

    %% STRATIFIED MEDIA SPECTRAL GREEN'S FUNCTION
    SGF(:, :, :, :, idx) = spectral_gf(dielectric.er(idx), wave.k0, ...
        kx, ky, v_tm, v_te, i_tm, i_te, 'E', 'M');
end

figure('Position', [250 250 750 400]);
for idx = 1 : 1 : length(dielectric.er)
    plot(kx / wave.k0, abs(SGF(:, :, 1, 2, idx)), ...
        'LineWidth', 2.0, 'DisplayName', ['\epsilon_{r} = ' ...
        num2str(dielectric.er(idx))]);
    hold on;
end
grid on;
xticks(0 : 0.25 : 1.5);
xlim([0 1.5]);
legend show;
legend('location', 'bestoutside');
xlabel('k_{x} / k_{0}');
ylabel('|G_{xy}|');
title(['Semi-Infinite Superstrate @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, h = ' num2str(bottom_medium.h * 1e3) ' mm, and z = h^{+}']);
saveas(gcf, 'figures\semi_infinite_superstrate.fig');
