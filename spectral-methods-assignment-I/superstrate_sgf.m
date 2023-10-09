close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
bottom_medium.er = 1;
bottom_medium.h = 15.6e-3;
dielectric.h  = 2.6e-3;
dielectric.er = 10;
top_medium.er = 1;
wave.f  = [8 8.5 9 9.5 10] * 1e9;
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
    
%% Z COORDINATE
z = ones(1, N) * (bottom_medium.h + dielectric.h);

ky = NaN(1, N, length(wave.f));
SGF = NaN(1, N, 3, 2, length(wave.f));
for idx = 1 : 1  : length(wave.f)
    %% WAVE VECTOR COMPONENTS
    ky(:, :, idx) = linspace(0, wave.k0(idx), N);
    kx = zeros( size(ky(:, :, idx), 1, 2) );
    krho = sqrt(kx .^ 2 + ky(:, :, idx) .^ 2);

    %% STRATIFIED MEDIA VOLTAGES AND CURRENTS
    [v_te, i_te, v_tm, i_tm] = stratified_media(wave.k0(idx), krho, z, ...
        'Superstrate', bottom_medium.h, dielectric.h, dielectric.er);

    %% STRATIFIED MEDIA SPECTRAL GREEN'S FUNCTION
    SGF(:, :, :, :, idx) = spectral_gf(top_medium.er, wave.k0(idx), ...
        kx, ky(:, :, idx), v_tm, v_te, i_tm, i_te, 'E', 'M');
end

figure('Position', [250 250 750 400]);
for idx = 1 : 1 : length(wave.f)
    plot(ky(:, :, idx) / wave.k0(idx), abs(SGF(:, :, 2, 1, idx)), ...
        'LineWidth', 2.0, 'DisplayName', ['f = ' ...
        num2str(wave.f(idx) * 1e-9) ' GHz']);
    hold on;
end
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('k_{x} / k_{0}');
ylabel('|G_{yx}|');
title(['Superstrate @ h = ' num2str(bottom_medium.h * 1e3) ...
    ' mm, h_{s} = ' num2str(dielectric.h * 1e3) 'mm, \epsilon_{r} = ' ...
    num2str(dielectric.er) ', and z = h + h_{s}^{+}']);
saveas(gcf, 'figures\superstrate.fig');
