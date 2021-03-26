import os
import numpy as np
import Core.Units
import Core.Constants
from decimal import Decimal
import matplotlib.pyplot as plt
import XSectParse.ParseXSectFile as x_parse
import XSectParse.PlotXSectFile as x_plot
from matplotlib.colors import LogNorm
import helper_functions as hf

# ==============================================================================
#                                  DECLARE DATA
# ==============================================================================

Li6_data = x_parse.parse('../data/Li6.txt')
B10_data = x_parse.parse('../data/B10.txt')
Cd113_data = x_parse.parse('../data/Cd113.txt')
Gd157_data = x_parse.parse('../data/Gd157.txt')

LiF_data = x_parse.parse('../data/LiF.txt')
LiF_99_perc_data = x_parse.parse('../data/LiF_99_perc.txt')
B4C_data = x_parse.parse('../data/B4C.txt')
B4C_99_perc_data = x_parse.parse('../data/B4C_99_perc.txt')
Cd_data = x_parse.parse('../data/Cd.txt')
Gd2O3_data = x_parse.parse('../data/Gd2O3.txt')

epoxy_gd_low_data = x_parse.parse('../data/epoxy_90_gd_10.txt')
epoxy_gd_high_data = x_parse.parse('../data/epoxy_10_gd_90.txt')

# Store elastic cross-sections in numpy arrays
Li6_array = np.transpose(np.array(Li6_data['procs']['hadElastic']))
B10_array = np.transpose(np.array(B10_data['procs']['hadElastic']))
Cd113_array = np.transpose(np.array(Cd113_data['procs']['hadElastic']))
Gd157_array = np.transpose(np.array(Gd157_data['procs']['hadElastic']))

LiF_array = np.transpose(np.array(LiF_data['procs']['hadElastic']))
LiF_99_perc_array = np.transpose(np.array(LiF_99_perc_data['procs']['hadElastic']))
B4C_array = np.transpose(np.array(B4C_data['procs']['hadElastic']))
B4C_99_perc_array = np.transpose(np.array(B4C_99_perc_data['procs']['hadElastic']))
Cd_array = np.transpose(np.array(Cd_data['procs']['hadElastic']))
Gd2O3_array = np.transpose(np.array(Gd2O3_data['procs']['hadElastic']))

epoxy_gd_low_array = np.transpose(np.array(epoxy_gd_low_data['procs']['hadElastic']))
epoxy_gd_high_array = np.transpose(np.array(epoxy_gd_high_data['procs']['hadElastic']))

# Extract densities
B4C_density = float(B4C_data['metadata']['NAtomsPerVolume [1/cm3]'])
B4C_99_perc_density = float(B4C_99_perc_data['metadata']['NAtomsPerVolume [1/cm3]'])
Cd_density = float(Cd_data['metadata']['NAtomsPerVolume [1/cm3]'])
Gd2O3_density = float(Gd2O3_data['metadata']['NAtomsPerVolume [1/cm3]'])

# Calculate densities
LiF_density = hf.get_number_density(2.635, 25.939, 2)
LiF_99_perc_density = hf.get_number_density(2.635, 25.939, 2)
epoxy_density_calculated = hf.get_number_density(1.18, 853.049, 123)
epoxy_gd_low_density = 0.1 * Gd2O3_density + 0.9 * epoxy_density_calculated
epoxy_gd_high_density = 0.9 * Gd2O3_density + 0.1 * epoxy_density_calculated

# ==============================================================================
#                           PLOT SCATTERING CROSS-SECTIONS
# ==============================================================================

# Scattering cross-sections (ISOTOPES)
hf.set_thick_labels(12)
fig = plt.figure()
plt.plot((Li6_array[0]/Core.Units.eV), Li6_array[1]/Core.Units.barn,
         color='blue', linestyle='-', label='$^{6}$Li')
plt.plot((B10_array[0]/Core.Units.eV), B10_array[1]/Core.Units.barn,
         color='red', linestyle='--', label='$^{10}$B')
plt.plot((Cd113_array[0]/Core.Units.eV), Cd113_array[1]/Core.Units.barn,
         color='green', linestyle='dotted', label='$^{113}$Cd')
plt.plot((Gd157_array[0]/Core.Units.eV), Gd157_array[1]/Core.Units.barn,
         color='black', linestyle='-.', label='$^{157}$Gd')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('$\sigma$ (barn)')
plt.title('Scattering cross-sections')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
ymin, ymax = 1e-2, 1e10

plt.fill_betweenx([ymin, ymax], 0.0001, 1, color='grey', alpha=0.2,
                 label=None, zorder=0)

plt.ylim(ymin, ymax)
#plt.legend(title='Isotopes')
fig.savefig('../output/isotope_scattering_vs_ev.pdf', bbox_inches='tight')

# Scattering cross-sections (MATERIALS)
fig = plt.figure()
hf.plot_cross_area(LiF_array, LiF_99_perc_array,
                'LiF ($^{6}$Li: 0.08-0.99)', 'blue', '-')

hf.plot_cross_area(B4C_array, B4C_99_perc_array,
                'B$_4$C ($^{10}$B: 0.20-0.99)', 'red', '--')

plt.plot((Cd_array[0]/Core.Units.eV), Cd_array[1]/Core.Units.barn,
         color='green', linestyle='dotted', label='Cd', zorder=4)

plt.plot((Gd2O3_array[0]/Core.Units.eV), Gd2O3_array[1]/Core.Units.barn,
         color='black', linestyle='-.', label='Gd$_2$O$_3$')

hf.plot_cross_area(epoxy_gd_low_array, epoxy_gd_high_array,
                'epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', 'orange', '-.')

plt.fill_betweenx([ymin, ymax], 0.0001, 1, color='grey', alpha=0.2,
                 label=None, zorder=0)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('$\sigma$ (barn)')
plt.title('Scattering cross-sections')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)

plt.ylim(ymin, ymax)
#plt.legend(title='Materials')
fig.savefig('../output/material_scattering_vs_ev.pdf', bbox_inches='tight')

# ==============================================================================
#                            PLOT CAPTURE PROBABILITIES
# ==============================================================================

# Declare arrays
arrays = [[LiF_array, LiF_99_perc_array],
          [B4C_array, B4C_99_perc_array],
          Cd_array, Gd2O3_array,
          [epoxy_gd_low_array, epoxy_gd_high_array]]
labels = ['LiF ($^{6}$Li: 0.08-0.99)', 'B$_4$C ($^{10}$B: 0.20-0.99)', 'Cd',
          'Gd$_2$O$_3$', 'epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)']
densities = [[LiF_density, LiF_99_perc_density],
             [B4C_density, B4C_99_perc_density],
             Cd_density, Gd2O3_density,
             [epoxy_gd_low_density, epoxy_gd_high_density]]
colors = ['blue', 'red', 'green', 'black', 'orange']
linestyles = ['-', '--', 'dotted', '-', '-.']
depth = 0.1
is_area_vec = [True, True, False, False, True]

# Plot - fixed thickness
fig = plt.figure()
for array_vec, label, density_vec, color, linestyle, is_area in zip(arrays, labels, densities, colors, linestyles, is_area_vec):
    if is_area:
        array_1, array_2 = array_vec
        d_1, d_2 = density_vec
        hf.plot_prob_area(array_1, array_2, label, color, linestyle, d_1, d_2, depth)
    else:
        array = array_vec
        density = density_vec
        plt.plot((array[0]/Core.Units.eV), 1 - np.exp(-(array[1]*(1e-2)*density*depth)),
                 color=color, linestyle=linestyle, label=label)
plt.title('Scattering probabilities - fixed thickness (1 mm)')
plt.xlabel('Energy (eV)')
plt.ylabel('Scattering probability')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Materials')
plt.xlim(1e-4, 1)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-3, 1)
fig.savefig('../output/scatering_probability_fixed_thickness_vs_ev.pdf', bbox_inches='tight')

# Plot - fixed energy
energy = 0.025
LiF_cross_idx = min(enumerate(LiF_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
LiF_99_cross_idx = min(enumerate(LiF_99_perc_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
B4C_cross_idx = min(enumerate(B4C_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
B4C_99_cross_idx = min(enumerate(B4C_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
Cd_cross_idx = min(enumerate(Cd_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
Gd2O3_cross_idx = min(enumerate(Gd2O3_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
epoxy_gd_low_idx = min(enumerate(epoxy_gd_low_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
epoxy_gd_high_idx = min(enumerate(epoxy_gd_high_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]

LiF_cross = LiF_array[1][LiF_cross_idx]
LiF_99_cross = LiF_99_perc_array[1][LiF_99_cross_idx]
B4C_cross = B4C_array[1][B4C_cross_idx]
B4C_99_cross = B4C_99_perc_array[1][B4C_99_cross_idx]
Cd_cross = Cd_array[1][Cd_cross_idx]
Gd2O3_cross = Gd2O3_array[1][Gd2O3_cross_idx]
epoxy_gd_low_cross = epoxy_gd_low_array[1][epoxy_gd_low_idx]
epoxy_gd_high_cross = epoxy_gd_high_array[1][epoxy_gd_high_idx]


cross_array = [[LiF_cross, LiF_99_cross],
               [B4C_cross, B4C_99_cross],
                Cd_cross, Gd2O3_cross,
               [epoxy_gd_low_cross, epoxy_gd_high_cross]]

thicknesses = np.arange(1e-4, 0.1, 0.00001)
fig = plt.figure()
plt.title('Scattering probability - fixed energy (0.025 eV)')
for cross_section_vec, density_vec, color, linestyle, label, is_area in zip(cross_array, densities, colors, linestyles, labels, is_area_vec):
    if is_area:
        cross_1, cross_2 = cross_section_vec
        d_1, d_2 = density_vec
        hf.plot_prob_area_thickness(cross_1, cross_2, label, color, linestyle, d_1, d_2, thicknesses)
    else:
        cross_section = cross_section_vec
        plt.plot(thicknesses*10, 1 - np.exp(-(cross_section*(1e-2)*density*thicknesses)),
                 color=color, linestyle=linestyle, label=label)
plt.legend(title='Material')
plt.xlabel('Thickness (mm)')
plt.ylabel('Scattering probability')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.001, 1)
#plt.ylim(0.00001, 1)
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
fig.savefig('../output/scattering_probability_fixed_energy_vs_depth.pdf', bbox_inches='tight')


# ==============================================================================
#                       PLOT THICKNESS VS ENERGY CURVES
# ==============================================================================

scattering_probability = 0.01
energy_start = 1e-5
energy_stop = 1
fig = plt.figure()
for array_vec, density_vec, color, linestyle, label, is_area in zip(arrays, densities, colors, linestyles, labels, is_area_vec):
    if is_area:
        array_1, array_2 = array_vec
        density_1, density_2 = density_vec
        hf.plot_thick_vs_energy_area(array_1, array_2, label, color, linestyle,
                                  density_1, density_2, energy_start, energy_stop,
                                  scattering_probability)
    else:
        array = array_vec
        density = density_vec
        # Extract cross-sections
        idx_start = min(enumerate(array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_start))[0]
        idx_stop = min(enumerate(array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_stop))[0]
        cross_sections = array[1][idx_start:idx_stop] * (1e-2)
        energies = array[0][idx_start:idx_stop]/Core.Units.eV
        # Extract thicknesses
        thicknesses = (- np.log(1 - scattering_probability) / (cross_sections * density)) * 10
        # Plot
        plt.plot(energies, thicknesses, color=color, linestyle=linestyle, label=label)
plt.title('Scattering probability = 1%')
plt.xlabel('Energy (eV)')
plt.ylabel('Thickness (mm)')
plt.xscale('log')
plt.yscale('log')
#plt.legend(title='Materials')
plt.xlim(0.00001, 1)
plt.ylim(1e-3, 1e5)
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
fig.savefig('../output/fixed_scattering_probability_energy_vs_thickness.pdf', bbox_inches='tight')



# ==============================================================================
#                                HELPER FUNCTION
# ==============================================================================

def meV_to_A(energy):
    return np.sqrt(81.81/energy)

def A_to_meV(wavelength):
    return (81.81 /(wavelength ** 2))
