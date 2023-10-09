close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
dielectric.h  = 4.5e-3;
dielectric.er = 6;
top_medium.er = 1;
wave.f  = [10 20] * 1e9;
N = 3001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
    
%% Z COORDINATE
z = ones(1, N) * dielectric.h;

kx = NaN(1, N, length(wave.f));
SGF = NaN(1, N, 3, 2, length(wave.f));
for idx = 1 : 1  : length(wave.f)
    %% WAVE VECTOR COMPONENTS
    kx(:, :, idx) = linspace(0, 3 * wave.k0(idx), N);
    ky = zeros( size(kx(:, :, idx), 1, 2) );
    krho = sqrt(kx(:, :, idx) .^ 2 + ky .^ 2);

    %% STRATIFIED MEDIA VOLTAGES AND CURRENTS
    [v_te, i_te, v_tm, i_tm] = stratified_media(wave.k0(idx), krho, z, ...
        'GroundSlab', dielectric.h, dielectric.er);

    %% STRATIFIED MEDIA SPECTRAL GREEN'S FUNCTION
    SGF(:, :, :, :, idx) = spectral_gf(top_medium.er, wave.k0(idx), ...
        kx(:, :, idx), ky, v_tm, v_te, i_tm, i_te, 'E', 'J');
end

figure('Position', [250 250 750 400]);
yyaxis left;
plot(kx(:, :, 1) / wave.k0(1), abs(SGF(:, :, 1, 1, 1)), ...
    'LineWidth', 2.0, 'DisplayName', ['f = ' ...
    num2str(wave.f(1) * 1e-9) ' GHz']);
hold on;
ylabel('|G_{xx}|');
yyaxis right;
plot(kx(:, :, 2) / wave.k0(2), abs(SGF(:, :, 1, 1, 2)), ...
    'LineWidth', 2.0, 'DisplayName', ['f = ' ...
    num2str(wave.f(2) * 1e-9) ' GHz']);
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('k_{x} / k_{0}');
ylabel('|G_{xx}|');
title(['Grounded Slab @ h = ' num2str(dielectric.h * 1e3) ...
    ' mm, \epsilon_{r} = ' num2str(dielectric.er) ', and z = h^{+}']);
saveas(gcf, 'figures\grounded_dielectric.fig');
