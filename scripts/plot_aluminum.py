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

Al_data = x_parse.parse('../data/Al.txt')
Al_array_scattering = np.transpose(np.array(Al_data['procs']['hadElastic']))
Al_array_absorption = hf.get_absorption_cross(Al_data, offset=0)
Al_density = float(Al_data['metadata']['NAtomsPerVolume [1/cm3]'])


# ==============================================================================
#                           PLOT CROSS-SECTIONS
# ==============================================================================

# Plot cross-sections
fig = plt.figure()
plt.plot((Al_array_scattering[0]/Core.Units.eV),
          Al_array_scattering[1]/Core.Units.barn,
          color='blue', linestyle='-', label='Scattering')
plt.plot((Al_array_absorption[0]/Core.Units.eV),
          Al_array_absorption[1]/Core.Units.barn,
          color='red', linestyle='--', label='Absorption')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('$\sigma$ (barn)')
plt.title('Cross-sections $^{27}$Al')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Process')
fig.savefig('../output/al_cross_sections_vs_ev.pdf', bbox_inches='tight')


# ==============================================================================
#                             PLOT PROBABILITIES
# ==============================================================================

# Declare arrays
depth_in_cm = 1.4
# Plot - fixed thickness
fig = plt.figure()
scat_prob = 1 - np.exp(-(Al_array_scattering[1]*(1e-2)*Al_density*depth_in_cm))
abs_prob = 1 - np.exp(-(Al_array_absorption[1]*(1e-2)*Al_density*depth_in_cm))
plt.plot((Al_array_scattering[0]/Core.Units.eV),
         scat_prob*100,
         color='blue', linestyle='-', label='Scattering')
plt.plot((Al_array_absorption[0]/Core.Units.eV),
         abs_prob*100,
         color='red', linestyle='--', label='Absorption')
plt.title('Aluminum - fixed thickness (14 mm)')
plt.xlabel('Energy (eV)')
plt.ylabel('Probability (%)')
plt.fill_betweenx([1e-5, 100], 0.00001, 1, color='grey', alpha=0.2,
                 label=None, zorder=0)
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Process')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-4, 1e2)
fig.savefig('../output/al_probability_fixed_thickness_vs_ev_%.2f_cm.pdf' % depth_in_cm, bbox_inches='tight')



# ==============================================================================
#                         PLOT PROBABILITIES - ZOOM +  WAVELENGTH
# ==============================================================================


# Declare arrays
depth = 1.25
# Fixed thickness
fig = plt.figure()
fig.set_figheight(5)
fig.set_figwidth(12)
wavelengths_scatter = hf.meV_to_A((Al_array_scattering[0]/Core.Units.eV)*1000)
probs_scatter = 1 - np.exp(-(Al_array_scattering[1]*(1e-2)*Al_density*depth))
plt.plot(wavelengths_scatter, probs_scatter*100,
         color='blue', linestyle='-', label='Scattering')
wavelengths_abs = hf.meV_to_A((Al_array_absorption[0]/Core.Units.eV)*1000)
probs_abs = 1 - np.exp(-(Al_array_absorption[1]*(1e-2)*Al_density*depth))
plt.plot(wavelengths_abs, probs_abs*100,
         color='red', linestyle='--', label='Absorption')
plt.title('Aluminum, 12.5 mm thickness (4 mm window + 8.5 mm substrates)')
plt.xlabel('Wavelength (Å)')
plt.ylabel('Probability (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Process')
plt.xlim(2, 20)
plt.ylim(0, 20)
fig.savefig('../output/al_probability_fixed_thickness_vs_wavelength.pdf', bbox_inches='tight')

# Save probablities in vectors and export (also fix correct number of elements)
np.savetxt('../output/al_scatter_prob_vs_wavelength.txt',
           np.transpose(np.array([wavelengths_scatter, probs_scatter])),
           delimiter=",",
           header='Wavelength (Å), Probability')
np.savetxt('../output/al_abs_prob_vs_wavelength.txt',
           np.transpose(np.array([wavelengths_abs, probs_abs])),
           delimiter=",",
           header='Wavelength (Å), Probability')


wavelengths = np.linspace(1, 7, 100)
attenuation = []
for wavelength in wavelengths:
    idx_abs = (np.abs(wavelengths_abs - wavelength)).argmin()
    idx_scatter = (np.abs(wavelengths_scatter - wavelength)).argmin()
    attenuation.append(probs_scatter[idx_scatter] * 0.5 + probs_abs[idx_abs])
    #print('Probability: %f' % (probs_scatter[idx_scatter] * 0.5 + probs_abs[idx_abs]))
attentuation = np.array(attenuation)
np.savetxt('../output/window_attenuation_vs_wavelength.txt',
           np.transpose(np.array([wavelengths, attenuation])),
           delimiter=",",
           header='Wavelength (Å), Attenuation (fraction)')



# ==============================================================================
#                                PLOT HEATMAP
# ==============================================================================

hf.set_thick_labels(12)
thicknesses = np.linspace(0, 15, 100)
energies = np.linspace(0.2, 200, 100) * 1e-3
heat_map = []
for i, energy in enumerate(energies):
    print('%d/%d' % (i+1, len(energies)))
    probabilities = np.empty(len(thicknesses), dtype='float')
    for j, thickness in enumerate(thicknesses):
        cross_idx = min(enumerate(Al_array_scattering[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
        cross_section = Al_array_scattering[1][cross_idx]
        prob_scatter = 1 - np.exp(-(cross_section*(1e-2)*Al_density*thickness*1e-1))
        probabilities[j] = prob_scatter * 100
    heat_map.append(probabilities)
heat_map = np.array(heat_map)
fig = plt.figure()
plt.imshow(heat_map, origin='lower', cmap='jet', aspect='auto',
           extent=[thicknesses[0], thicknesses[-1], energies[0], energies[-1]])
cbar = plt.colorbar()
cbar.set_label('Probability (%)')
plt.xlabel('Thickness (mm)')
plt.ylabel('Energy (eV)')
plt.title('Aluminum - scattering probability')
plt.scatter([8.5], [0.020], marker='x', color='black', label='CSPEC', s=[50])
plt.scatter([10.5], [0.168], marker='+', color='black', label='T-REX', s=[60])
#plt.plot([8.5, 8.5], [2e-4, 0.020], linestyle='-', color='black', label='CSPEC')
#plt.plot([10.5, 10.5], [0.002, 0.168], linestyle='--', color='black', label='T-REX')
plt.legend(title='Detector', loc=2)
fig.savefig('../output/al_scatter_function_of_thickness_and_energy.pdf', bbox_inches='tight')
