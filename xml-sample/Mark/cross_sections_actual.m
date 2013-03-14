% generate macroscopic XS's and other input data from Mark's XS data

clear all
close all
clc

load '/home/will/Documents/NSE_Research/openmoc/xml-sample/Mark/SFRgroupXSmod.mat';

% Input a value for F
F = 1;

% material densities (atoms / barn / cm)
U235_density_assembly = .0029946 / F^2;
U238_density_assembly = .0269514 / F^2;
NA23_density_assembly = .009234 * (1 / F^2) * (F^2 - 0.62) / 0.38;
NA23_density_blanket  = .0243;

% macroscopic cross sections material 1 (abs, total, fis, nu_fis, chi, scat)
% assembly
xs_fis_assembly = (U235_density_assembly * U235f + U238_density_assembly * U238f)';
xs_abs_assembly = (U235_density_assembly * U235a + U238_density_assembly * U238a...
    + NA23_density_assembly * NA23a)';% + xs_fis_assembly; %SHOULD THIS BE HERE?
xs_nu_fis_assembly = nu' .* xs_fis_assembly;
xs_scat_assembly = (U235_density_assembly * U235s + U238_density_assembly * U238s...
    + NA23_density_assembly * NA23s);
chi_assembly = chi';
xs_tot_assembly = xs_abs_assembly + sum(xs_scat_assembly,2)';

% blanket
xs_abs_blanket = (NA23_density_blanket * NA23a)';
xs_fis_blanket = zeros(1,33);
xs_nu_fis_blanket = zeros(1,33);
xs_scat_blanket = (NA23_density_blanket * NA23s);
chi_blanket = zeros(1,33);
xs_tot_blanket = xs_abs_blanket + sum(xs_scat_blanket,2)';
