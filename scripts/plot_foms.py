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

# Import data
LiF_data = x_parse.parse('../data/LiF.txt')
LiF_99_perc_data = x_parse.parse('../data/LiF_99_perc.txt')
B4C_data = x_parse.parse('../data/B4C.txt')
B4C_97_perc_data = x_parse.parse('../data/B4C_97_perc.txt')
B4C_99_perc_data = x_parse.parse('../data/B4C_99_perc.txt')
Cd_data = x_parse.parse('../data/Cd.txt')
Gd2O3_data = x_parse.parse('../data/Gd2O3.txt')
epoxy_gd_low_data = x_parse.parse('../data/epoxy_90_gd_10.txt')
epoxy_gd_medium_data = x_parse.parse('../data/epoxy_50_gd_50.txt')
epoxy_gd_high_data = x_parse.parse('../data/epoxy_10_gd_90.txt')
Al_data = x_parse.parse('../data/Al.txt')
boral_data = x_parse.parse('../data/boral.txt')

# Extract densities
B4C_density = float(B4C_data['metadata']['NAtomsPerVolume [1/cm3]'])
B4C_97_perc_density = float(B4C_97_perc_data['metadata']['NAtomsPerVolume [1/cm3]'])
B4C_99_perc_density = float(B4C_99_perc_data['metadata']['NAtomsPerVolume [1/cm3]'])
Cd_density = float(Cd_data['metadata']['NAtomsPerVolume [1/cm3]'])
Gd2O3_density = float(Gd2O3_data['metadata']['NAtomsPerVolume [1/cm3]'])
LiF_density = hf.get_number_density(2.635, 25.939, 2)
LiF_99_perc_density = hf.get_number_density(2.635, 25.939, 2)
epoxy_density_calculated = hf.get_number_density(1.18, 853.049, 123)
epoxy_gd_low_density = 0.9 * Gd2O3_density + 0.1 * epoxy_density_calculated
epoxy_gd_medium_density = 0.5 * Gd2O3_density + 0.5 * epoxy_density_calculated
epoxy_gd_high_density = 0.9 * Gd2O3_density + 0.1 * epoxy_density_calculated
Al_density = float(Al_data['metadata']['NAtomsPerVolume [1/cm3]'])
boral_density = 0.67 * B4C_density + 0.33 * Al_density

# Put everything together in one array
data_array = np.array([['LiF_1', LiF_data, LiF_density],
                      ['LiF_2', LiF_99_perc_data, LiF_99_perc_density],
                      ['B4C_1', B4C_data, B4C_density],
                      ['B4C_2', B4C_99_perc_data, B4C_99_perc_density],
                      ['Cd', Cd_data, Cd_density],
                      ['Gd2O3', Gd2O3_data, Gd2O3_density],
                      ['epoxyGd2O3_1', epoxy_gd_low_data, epoxy_gd_low_density],
                      ['epoxyGd2O3_2', epoxy_gd_high_data, epoxy_gd_high_density],
                      ['boral', boral_data, boral_density]
                      ])

title_names = ['LiF ($^6$Li: 0.08)', 'LiF ($^6$Li: 0.99)',
               'B$_4$C ($^{10}$B: 0.20)', 'B$_4$C ($^{10}$B: 0.99)',
               'Cd', 'Gd$_2$O$_3$', 'epoxy-Gd$_2$O$_3$ (w/w ratio: 0.1-0.9)',
               'epoxy-Gd$_2$O$_3$ (w/w ratio: 0.9-0.1)', 'Boral']

# ==============================================================================
#                    GET ALL HEAT MAPS
# ==============================================================================

def get_all_heat_maps(data_array):
    @plt.FuncFormatter
    def fake_log(x, pos):
        'The two args are the value and tick position'
        return r'$10^{%d}$' % (x)
    # Declare parameters
    hf.set_thick_labels(15)
    number_elements = 400
    thicknesses_in_cm = np.logspace(-3, 2, number_elements) * 0.1  # in cm
    energies_in_eV = np.logspace(-5, 0, number_elements)
    for i, (material, data, density) in enumerate(data_array):
        absorption_array = hf.get_absorption_cross(data)
        scattering_array = hf.get_scattering_cross(data)
        heat_map = []
        print(material)
        for j, energy in enumerate(energies_in_eV):
            print('%d/%d' % (j+1, len(energies_in_eV)))
            probabilities = np.empty(len(thicknesses_in_cm), dtype='float')
            for k, thickness in enumerate(thicknesses_in_cm):
                # Extract cross-sections
                cross_idx = min(enumerate(absorption_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy))[0]
                cross_section_abs = absorption_array[1][cross_idx]
                cross_section_scat = scattering_array[1][cross_idx]
                # Get scattering and absorption probabilities
                p_absorption = 1 - np.exp(-(cross_section_abs*(1e-2)*density*thickness))
                p_scattering = 1 - np.exp(-(cross_section_scat*(1e-2)*density*thickness))
                # Get transmission
                p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
                probabilities[k] = p_transmission * 100 + 1e-20 # a very small number is added just to keep the probability from being 0
            heat_map.append(probabilities)
        heat_map = np.transpose(np.array(heat_map))
        fig, (ax1) = plt.subplots()
        plt.imshow(heat_map, origin='lower', cmap='jet', aspect='auto', norm=LogNorm(),
                   extent=[-5, 0, -3, 2], vmin=1e-4, vmax=1e2)
        ax1.imshow(heat_map, origin='lower', cmap='jet', aspect='auto', norm=LogNorm(),
                   extent=[-5, 0, -3, 2], vmin=1e-4, vmax=1e2)
        cbar = plt.colorbar(extend='min')
        cbar.set_label('Transmission (%)')
        ax1.set_xlabel('Energy (eV)')
        ax1.set_ylabel('Thickness (mm)')
        ax1.xaxis.set_major_formatter(fake_log)
        ax1.yaxis.set_major_formatter(fake_log)
        ax1.grid(True, which='major', linestyle='--', zorder=5)
        #plt.xscale('log')
        #plt.yscale('log')
        ax1.set_title(title_names[i])
        fig.set_figheight(6)
        fig.set_figwidth(10)
        fig.savefig('../output/trans_function_of_thickness_and_energy_%s.pdf' % material, bbox_inches='tight')
    print('Done!')

#get_all_heat_maps(data_array)
#print('Jaadå')

# ==============================================================================
#                          GET ABSORPTION AND SCATTERING
# ==============================================================================

def get_absorption_and_scattering(data_array):
    # Declare parameters
    data_dict = {'LiF_1': None,
                 'LiF_2': None,
                 'B4C_1': None,
                 'B4C_2': None,
                 'Cd': None,
                 'Gd2O3': None,
                 'epoxyGd2O3_1': None,
                 'epoxyGd2O3_2': None,
                 'boral': None}
    energies = np.logspace(-5, 0, 1000)
    energy_start, energy_stop = 1e-5, 1
    capture_probability = 0.99
    for material, data, density in data_array:
        scattering_array = hf.get_scattering_cross(data)
        absorption_array = hf.get_absorption_cross(data)
        # Extract cross-sections
        idx_start = min(enumerate(absorption_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_start))[0]
        idx_stop = min(enumerate(absorption_array[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_stop))[0]
        cross_sections = absorption_array[1][idx_start:idx_stop] * (1e-2)
        energies = absorption_array[0][idx_start:idx_stop]/Core.Units.eV
        # Extract thicknesses
        thicknesses = (-np.log(1 - capture_probability) / (cross_sections * density)) # in cm
        # Extract corresponding scattering probabilities
        scattering_probs = 1 - np.exp(-(scattering_array[1][idx_start:idx_stop]*(1e-2)*density*thicknesses))
        data_dict[material] = [energies, thicknesses*10, scattering_probs*100]

    # Plot data
    thicknesses_in_mm = np.logspace(-2, 2, 1000)[1:]
    hf.set_thick_labels(12)
    fig = plt.figure()
    fig.suptitle('Absorption = 99 %')
    plt.subplot(1, 2, 1)
    plt.plot(data_dict['Cd'][0], data_dict['Cd'][1], color='green', linestyle='dotted', label='Cd')
    plt.plot(data_dict['Gd2O3'][0], data_dict['Gd2O3'][1], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(data_dict['LiF_1'][0], data_dict['LiF_1'][1], color='blue', linestyle='solid', label=None)
    plt.plot(data_dict['LiF_2'][0], data_dict['LiF_2'][1], color='blue', linestyle='solid', label=None)
    plt.fill_between(data_dict['LiF_2'][0], data_dict['LiF_2'][1], data_dict['LiF_1'][1], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(data_dict['B4C_1'][0], data_dict['B4C_1'][1], color='red', linestyle='--', label=None)
    plt.plot(data_dict['B4C_2'][0], data_dict['B4C_2'][1], color='red', linestyle='--', label=None)
    plt.fill_between(data_dict['B4C_2'][0], data_dict['B4C_2'][1], data_dict['B4C_1'][1], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(data_dict['epoxyGd2O3_1'][0], data_dict['epoxyGd2O3_1'][1], color='orange', linestyle='-.', label=None)
    plt.plot(data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_2'][1], color='orange', linestyle='-.', label=None)
    plt.fill_between(data_dict['epoxyGd2O3_1'][0], data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)
    plt.title('Thickness')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Thickness (mm)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(title='Materials')
    plt.xlim(0.0001, 1)
    plt.ylim(1e-3, 1e5)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)

    plt.subplot(1, 2, 2)
    plt.plot(data_dict['Cd'][0], data_dict['Cd'][2], color='green', linestyle='dotted', label='Cd')
    plt.plot(data_dict['Gd2O3'][0], data_dict['Gd2O3'][2], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(data_dict['LiF_1'][0], data_dict['LiF_1'][2], color='blue', linestyle='solid', label=None)
    plt.plot(data_dict['LiF_2'][0], data_dict['LiF_2'][2], color='blue', linestyle='solid', label=None)
    plt.fill_between(data_dict['LiF_2'][0], data_dict['LiF_2'][2], data_dict['LiF_1'][2], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(data_dict['B4C_1'][0], data_dict['B4C_1'][2], color='red', linestyle='--', label=None)
    plt.plot(data_dict['B4C_2'][0], data_dict['B4C_2'][2], color='red', linestyle='--', label=None)
    plt.fill_between(data_dict['B4C_2'][0], data_dict['B4C_2'][2], data_dict['B4C_1'][2], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(data_dict['epoxyGd2O3_1'][0], data_dict['epoxyGd2O3_1'][2], color='orange', linestyle='-.', label=None)
    plt.plot(data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_2'][2], color='orange', linestyle='-.', label=None)
    plt.fill_between(data_dict['epoxyGd2O3_1'][0], data_dict['epoxyGd2O3_2'][2], data_dict['epoxyGd2O3_1'][2], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)
    plt.title('Scattering')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Scattering probability (%)')
    plt.xscale('log')
    plt.yscale('log')
    #plt.legend(title='Materials')
    plt.xlim(0.0001, 1)
    plt.ylim(1e-1, 1e2)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)


    fig.set_figheight(5)
    fig.set_figwidth(12)
    plt.tight_layout()
    fig.savefig('../output/absorption_99_and_scattering.pdf', bbox_inches='tight')
    #print('Thicknesses Cd (absorption):')
    #print(data_dict['Cd'][1])

#get_absorption_and_scattering(data_array)
#print('Vi är klara!!')


# ==============================================================================
#       GET ALBEDO SATURATION LEVEL AND THICKNESS AS A FUNCTION OF ENERGY
# ==============================================================================

def get_sat_level_and_thickness(data_array):
    # Declare parameters
    number_elements = 500
    thicknesses_in_cm = np.logspace(-6, 2, number_elements) * 0.1  # in cm
    energies_in_eV = np.logspace(-5, 0, number_elements)
    data_dict = {'LiF_1': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'LiF_2': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'B4C_1': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'B4C_2': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'Cd':    np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'Gd2O3': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'epoxyGd2O3_1': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'epoxyGd2O3_2': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)]),
                 'boral': np.array([np.empty(number_elements), np.empty(number_elements), np.empty(number_elements)])
                 }

    for i, energy in enumerate(energies_in_eV):
        print(energy)
        for j, (material, data, density) in enumerate(data_array):
            scattering_array = hf.get_scattering_cross(data)
            absorption_array = hf.get_absorption_cross(data)
            # Extract cross-sections
            sigma_scatt = hf.get_sigma(scattering_array, energy)
            sigma_abs = hf.get_sigma(absorption_array, energy)
            # Get scattering and absorption probabilities
            p_scattering = hf.get_probability(sigma_scatt, density, thicknesses_in_cm)
            p_absorption = hf.get_probability(sigma_abs, density, thicknesses_in_cm)
            # Get transmission and albedo probabilities
            p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
            p_albedo_elements = ((1 - p_absorption[1:]) ** 2) * 0.5 * np.diff(p_scattering)
            p_albedo = []
            for k in range(0, len(p_albedo_elements)):
                p_albedo.append(sum(p_albedo_elements[:k+1]))
            p_albedo = np.array(p_albedo)
            # Get the first index of where the albedo level saturates
            albedo_sat_level = max(p_albedo) * 100
            albedo_0_95_sat_idx = np.where(p_albedo*100 >= albedo_sat_level*0.95)[0][0]
            albedo_sat_thick = thicknesses_in_cm[albedo_0_95_sat_idx]
            # Store in data_dict
            data_dict[material][0][i] = albedo_sat_thick * 10 # in mm
            data_dict[material][1][i] = albedo_sat_level * 0.95
            data_dict[material][2][i] = (p_transmission[albedo_0_95_sat_idx]*100)/(albedo_sat_level * 0.95)
    # Plot data
    hf.set_thick_labels(12)
    fig = plt.figure()
    plt.plot(energies_in_eV, data_dict['Cd'][0], color='green', linestyle='dotted', label='Cd')
    plt.plot(energies_in_eV, data_dict['Gd2O3'][0], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(energies_in_eV, data_dict['LiF_1'][0], color='blue', linestyle='solid', label=None)
    plt.plot(energies_in_eV, data_dict['LiF_2'][0], color='blue', linestyle='solid', label=None)
    plt.fill_between(energies_in_eV, data_dict['LiF_2'][0], data_dict['LiF_1'][0], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['B4C_1'][0], color='red', linestyle='--', label=None)
    plt.plot(energies_in_eV, data_dict['B4C_2'][0], color='red', linestyle='--', label=None)
    plt.fill_between(energies_in_eV, data_dict['B4C_2'][0], data_dict['B4C_1'][0], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_1'][0], color='orange', linestyle='-.', label=None)
    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_2'][0], color='orange', linestyle='-.', label=None)
    plt.fill_between(energies_in_eV, data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_1'][0], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

    plt.xlabel('Energy (eV)')
    plt.ylabel('Thickness (mm)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Thickness')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlim(1e-4, 1)
    plt.ylim(1e-3, 1e2)
    plt.legend(title='Material', loc=2)
    fig.savefig('../output/albedo_sat_thick_vs_energy.pdf' % energy, bbox_inches='tight')

    fig = plt.figure()
    plt.plot(energies_in_eV, data_dict['Cd'][1], color='green', linestyle='dotted', label='Cd')
    plt.plot(energies_in_eV, data_dict['Gd2O3'][1], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(energies_in_eV, data_dict['LiF_1'][1], color='blue', linestyle='solid', label=None)
    plt.plot(energies_in_eV, data_dict['LiF_2'][1], color='blue', linestyle='solid', label=None)
    plt.fill_between(energies_in_eV, data_dict['LiF_2'][1], data_dict['LiF_1'][1], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['B4C_1'][1], color='red', linestyle='--', label=None)
    plt.plot(energies_in_eV, data_dict['B4C_2'][1], color='red', linestyle='--', label=None)
    plt.fill_between(energies_in_eV, data_dict['B4C_2'][1], data_dict['B4C_1'][1], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_1'][1], color='orange', linestyle='-.', label=None)
    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_2'][1], color='orange', linestyle='-.', label=None)
    plt.fill_between(energies_in_eV, data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

    plt.xlabel('Energy (eV)')
    plt.ylabel('Albedo (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Level')
    plt.xlim(1e-4, 1)
    plt.ylim(1e-4, 1e2)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    fig.savefig('../output/albedo_sat_level_vs_energy.pdf' % energy, bbox_inches='tight')

    fig = plt.figure()
    plt.plot(energies_in_eV, data_dict['Cd'][2], color='green', linestyle='dotted', label='Cd')
    plt.plot(energies_in_eV, data_dict['Gd2O3'][2], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(energies_in_eV, data_dict['LiF_1'][2], color='blue', linestyle='solid', label=None)
    plt.plot(energies_in_eV, data_dict['LiF_2'][2], color='blue', linestyle='solid', label=None)
    plt.fill_between(energies_in_eV, data_dict['LiF_2'][2], data_dict['LiF_1'][2], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['B4C_1'][2], color='red', linestyle='--', label=None)
    plt.plot(energies_in_eV, data_dict['B4C_2'][2], color='red', linestyle='--', label=None)
    plt.fill_between(energies_in_eV, data_dict['B4C_2'][2], data_dict['B4C_1'][2], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_1'][2], color='orange', linestyle='-.', label=None)
    plt.plot(energies_in_eV, data_dict['epoxyGd2O3_2'][2], color='orange', linestyle='-.', label=None)
    plt.fill_between(energies_in_eV, data_dict['epoxyGd2O3_2'][2], data_dict['epoxyGd2O3_1'][2], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

    plt.xlabel('Energy (eV)')
    plt.ylabel('Transmission/Albedo')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Transmission/Albedo at saturation thickness')
    plt.xlim(1e-4, 1)
    #plt.ylim(1e-4, 1e2)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    fig.savefig('../output/albedo_sat_trans_over_albedo_vs_energy.pdf' % energy, bbox_inches='tight')
    # Save albedo saturation levels in text
    for material, (__, albedo_sat, __) in data_dict.items():
        np.savetxt('../output/%s_albedo_saturation_level.txt' % material,
                   np.transpose(np.array([energies_in_eV, albedo_sat])),
                   delimiter=",",
                   header='Energies (eV), Albedo saturation level (%)')




#get_sat_level_and_thickness(data_array)
#print('Done with new stuff!')

# ==============================================================================
#                  GET TRANSMISSION AND ALBEDO FOR A FIXED ENERGY
# ==============================================================================

def get_transmission_and_albedo(data_array, energy, title):
    thicknesses_in_cm = np.logspace(-7, 2, 1000) * 0.1  # in cm
    data_dict = {'LiF_1': None,
                 'LiF_2': None,
                 'B4C_1': None,
                 'B4C_2': None,
                 'Cd': None,
                 'Gd2O3': None,
                 'epoxyGd2O3_1': None,
                 'epoxyGd2O3_2': None,
                 'boral': None}

    for material, data, density in data_array:
        scattering_array = hf.get_scattering_cross(data)
        absorption_array = hf.get_absorption_cross(data)
        # Extract cross-sections
        sigma_scatt = hf.get_sigma(scattering_array, energy)
        sigma_abs = hf.get_sigma(absorption_array, energy)
        # Get scattering and absorption probabilities
        p_scattering = hf.get_probability(sigma_scatt, density, thicknesses_in_cm)
        p_absorption = hf.get_probability(sigma_abs, density, thicknesses_in_cm)
        # Get transmission and albedo probabilities
        p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
        p_albedo_elements = ((1 - p_absorption[1:]) ** 2) * 0.5 * np.diff(p_scattering)
        p_albedo = []
        for j in range(0, len(p_albedo_elements)):
            p_albedo.append(sum(p_albedo_elements[:j+1]))
        p_albedo = np.array(p_albedo)
        # Get the first index of where the albedo level saturates
        #print(material)
        #print(np.diff(p_albedo*100))
        #albedo_sat_idx = np.where(np.diff(p_albedo*100) <= 1e-4)[0][0]
        albedo_sat_level = max(p_albedo) * 100
        albedo_0_95_sat_idx = np.where(p_albedo*100 >= albedo_sat_level*0.95)[0][0]
        # Store in data_dict
        data_dict[material] = [p_transmission[1:]*100, p_albedo*100, albedo_0_95_sat_idx]

    # Plot data
    thicknesses_in_mm = np.logspace(-7, 2, 1000)[1:]
    hf.set_thick_labels(12)
    fig = plt.figure()
    fig.suptitle(title)
    plt.subplot(1, 2, 1)
    plt.plot(thicknesses_in_mm, data_dict['Cd'][0], color='green', linestyle='dotted', label='Cd')
    plt.plot(thicknesses_in_mm, data_dict['Gd2O3'][0], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(thicknesses_in_mm, data_dict['LiF_1'][0], color='blue', linestyle='solid', label=None)
    plt.plot(thicknesses_in_mm, data_dict['LiF_2'][0], color='blue', linestyle='solid', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['LiF_2'][0], data_dict['LiF_1'][0], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['B4C_1'][0], color='red', linestyle='--', label=None)
    plt.plot(thicknesses_in_mm, data_dict['B4C_2'][0], color='red', linestyle='--', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['B4C_2'][0], data_dict['B4C_1'][0], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][0], color='orange', linestyle='-.', label=None)
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], color='orange', linestyle='-.', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_1'][0], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

    #plt.plot(thicknesses_in_mm, data_dict['boral'][0], color='grey', linestyle='-', label='Boral')

    plt.xlabel('Thickness (mm)')
    plt.ylabel('Transmission (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Transmission')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    plt.legend(title='Material', loc=3)

    plt.subplot(1, 2, 2)
    plt.plot(thicknesses_in_mm, data_dict['Cd'][1], color='green', linestyle='dotted', label='Cd')
    #plt.axvline(thicknesses_in_mm[data_dict['Cd'][2]-1], color='green', label=None)
    plt.plot(thicknesses_in_mm, data_dict['Gd2O3'][1], color='black', linestyle='solid', label='Gd$_2$O$_3$')
    #plt.axvline(thicknesses_in_mm[data_dict['Gd2O3'][2]-1], color='black', label=None)

    plt.plot(thicknesses_in_mm, data_dict['LiF_1'][1], color='blue', linestyle='solid', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['LiF_1'][2]-1], color='blue', label=None)
    plt.plot(thicknesses_in_mm, data_dict['LiF_2'][1], color='blue', linestyle='solid', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['LiF_2'][2]-1], color='blue', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['LiF_2'][1], data_dict['LiF_1'][1], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['B4C_1'][1], color='red', linestyle='--', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['B4C_1'][2]-1], color='red', label=None)
    plt.plot(thicknesses_in_mm, data_dict['B4C_2'][1], color='red', linestyle='--', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['B4C_2'][2]-1], color='red', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['B4C_2'][1], data_dict['B4C_1'][1], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][1], color='orange', linestyle='-.', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['epoxyGd2O3_1'][2]-1], color='orange', label=None)
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], color='orange', linestyle='-.', label=None)
    #plt.axvline(thicknesses_in_mm[data_dict['epoxyGd2O3_2'][2]-1], color='orange', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

    #plt.plot(thicknesses_in_mm, data_dict['boral'][1], color='grey', linestyle='-', label='Boral')

    plt.xlabel('Thickness (mm)')
    plt.ylabel('Albedo (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Albedo')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    fig.set_figheight(5)
    fig.set_figwidth(12)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    plt.tight_layout()
    fig.savefig('../output/albedo_and_transmission_at_%.4f_eV.pdf' % energy, bbox_inches='tight')

    # Save transmission and albedo in text
    for material, (p_trans, p_albedo, __) in data_dict.items():
        np.savetxt('../output/%s_%.2f_meV.txt' % (material, energy*1000),
                   np.transpose(np.array([thicknesses_in_mm, p_trans, p_albedo])),
                   delimiter=",",
                   header='Thicknesses (mm), Transmission (%), Albedo (%)')


def get_transmission_and_albedo_epoxy(data_array, energy, title):
    thicknesses_in_cm = np.logspace(-2, 2, 1000) * 0.1  # in cm
    data_dict = {'LiF_1': None,
                 'LiF_2': None,
                 'B4C_1': None,
                 'B4C_2': None,
                 'Cd': None,
                 'Gd2O3': None,
                 'epoxyGd2O3_1': None,
                 'epoxyGd2O3_2': None,
                 'boral': None}

    for material, data, density in data_array:
        scattering_array = hf.get_scattering_cross(data)
        absorption_array = hf.get_absorption_cross(data)
        # Extract cross-sections
        sigma_scatt = hf.get_sigma(scattering_array, energy)
        sigma_abs = hf.get_sigma(absorption_array, energy)
        # Get scattering and absorption probabilities
        p_scattering = hf.get_probability(sigma_scatt, density, thicknesses_in_cm)
        p_absorption = hf.get_probability(sigma_abs, density, thicknesses_in_cm)
        # Get transmission and albedo probabilities
        p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
        p_albedo_elements = ((1 - p_absorption[1:]) ** 2) * 0.5 * np.diff(p_scattering)
        p_albedo = []
        for j in range(0, len(p_albedo_elements)):
            p_albedo.append(sum(p_albedo_elements[:j+1]))
        p_albedo = np.array(p_albedo)
        # Store in data_dict
        data_dict[material] = [p_transmission[1:]*100, p_albedo*100]

    # Plot data
    thicknesses_in_mm = np.logspace(-2, 2, 1000)[1:]
    hf.set_thick_labels(12)
    fig = plt.figure()
    fig.suptitle(title)
    plt.subplot(1, 2, 1)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][0], color='red', linestyle='-', label='epoxy-Gd$_2$O$_3$, 0.9-0.1 (w/w)')
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], color='blue', linestyle='dashed', label='epoxy-Gd$_2$O$_3$, 0.1-0.9 (w/w)')
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_1'][0], color='grey',
                     label=None, alpha=0.3)

    plt.xlabel('Thickness (mm)')
    plt.ylabel('Transmission (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Transmission')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    plt.legend(title='Material', loc=3)

    plt.subplot(1, 2, 2)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][1], color='red', linestyle='-', label='epoxy-Gd$_2$O$_3$, 0.9-0.1 (w/w)')
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], color='blue', linestyle='dashed', label='epoxy-Gd$_2$O$_3$, 0.1-0.9 (w/w)')
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='grey',
                     label=None, alpha=0.3)

    plt.xlabel('Thickness (mm)')
    plt.ylabel('Albedo (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Albedo')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    fig.set_figheight(5)
    fig.set_figwidth(12)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    plt.tight_layout()
    fig.savefig('../output/albedo_and_transmission_at_%.6f_eV.pdf' % energy, bbox_inches='tight')

get_transmission_and_albedo(data_array, 0.0001, 'Energy = 0.1 meV (28.6 Å)')
#get_transmission_and_albedo(data_array, 0.16696, 'Energy = 167 meV (0.7 Å)')
#get_transmission_and_albedo(data_array, 0.00021, 'Energy = 0.2 meV (20.0 Å)')
#get_transmission_and_albedo(data_array, 0.02525, 'Energy = 25 meV (1.8 Å)')
#get_transmission_and_albedo(data_array, 0.16696, 'Energy = 167 meV (0.7 Å)')
#get_transmission_and_albedo(data_array, 0.00001, 'Energy = 0.01 meV (90.45 Å)')
#get_transmission_and_albedo(data_array, 1.017, 'Energy = 1017 meV (0.286 Å)')
print('Wooop woop it has been done')
#get_transmission_and_albedo(data_array, 0.00001, 'Energy = 1e-5 eV')
#get_transmission_and_albedo(data_array, 0.0001, 'Energy = 1e-4 eV')
#get_transmission_and_albedo(data_array, 0.001, 'Energy = 1e-3 eV')
#get_transmission_and_albedo(data_array, 0.01, 'Energy = 1e-2 eV')
#get_transmission_and_albedo(data_array, 0.1, 'Energy = 1e-1 eV')
#get_transmission_and_albedo(data_array, 1, 'Energy = 1 eV')


#get_transmission_and_albedo_epoxy(data_array, 0.020, 'Energy = 20 meV (2 Å)')
#get_transmission_and_albedo_epoxy(data_array, 0.0002, 'Energy = 0.2 meV (20 Å)')
#get_transmission_and_albedo_epoxy(data_array, 0.0008, 'Energy = 0.8 meV (10 Å)')

# ==============================================================================
#                                GET FIGURES-OF-MERIT
# ==============================================================================

def get_foms(data_array, energy, p_backscattering=0.5):
    thicknesses_in_cm = np.logspace(-7, 2, 1000) * 0.1  # in cm
    data_dict = {'LiF_1': None,
                 'LiF_2': None,
                 'B4C_1': None,
                 'B4C_2': None,
                 'Cd': None,
                 'Gd2O3': None,
                 'epoxyGd2O3_1': None,
                 'epoxyGd2O3_2': None}

    for material, data, density in data_array:
        scattering_array = hf.get_scattering_cross(data)
        absorption_array = hf.get_absorption_cross(data)
        # Extract cross-sections
        sigma_scatt = hf.get_sigma(scattering_array, energy)
        sigma_abs = hf.get_sigma(absorption_array, energy)
        # Get scattering and absorption probabilities
        p_scattering = hf.get_probability(sigma_scatt, density, thicknesses_in_cm)
        p_absorption = hf.get_probability(sigma_abs, density, thicknesses_in_cm)
        # Get transmission and albedo probabilities
        p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
        p_albedo_elements = ((1 - p_absorption[1:]) ** 2) * 0.5 * np.diff(p_scattering)
        p_albedo = []
        for j in range(0, len(p_albedo_elements)):
            p_albedo.append(sum(p_albedo_elements[:j+1]))
        p_albedo = np.array(p_albedo)
        # Get figures-of-merit
        p_side = p_transmission[1:]
        p_back = p_albedo + (p_transmission[1:] ** 2) * p_backscattering
        p_internal = p_transmission[1:] + p_albedo
        # Store in data_dict
        data_dict[material] = [p_side*100, p_back*100, p_internal*100]

    # Plot data
    thicknesses_in_mm = np.logspace(-7, 2, 1000)[1:]
    hf.set_thick_labels(12)

    # SIDE
    fig = plt.figure()
    plt.plot(thicknesses_in_mm, data_dict['Cd'][0], color='green', linestyle='dotted', label='Cd')
    plt.plot(thicknesses_in_mm, data_dict['Gd2O3'][0], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(thicknesses_in_mm, data_dict['LiF_1'][0], color='blue', linestyle='solid', label=None)
    plt.plot(thicknesses_in_mm, data_dict['LiF_2'][0], color='blue', linestyle='solid', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['LiF_2'][0], data_dict['LiF_1'][0], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['B4C_1'][0], color='red', linestyle='--', label=None)
    plt.plot(thicknesses_in_mm, data_dict['B4C_2'][0], color='red', linestyle='--', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['B4C_2'][0], data_dict['B4C_1'][0], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][0], color='orange', linestyle='-.', label=None)
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], color='orange', linestyle='-.', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_1'][0], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)
    plt.xlabel('Thickness (mm)')
    plt.ylabel('$P_{transmission}$ (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('$fom_{side}$')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.legend(title='Material')
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    fig.savefig('../output/fom_side_at_%.4f_eV.pdf' % energy, bbox_inches='tight')

    # BACK
    fig = plt.figure()
    plt.plot(thicknesses_in_mm, data_dict['Cd'][1], color='green', linestyle='dotted', label='Cd')
    plt.plot(thicknesses_in_mm, data_dict['Gd2O3'][1], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(thicknesses_in_mm, data_dict['LiF_1'][1], color='blue', linestyle='solid', label=None)
    plt.plot(thicknesses_in_mm, data_dict['LiF_2'][1], color='blue', linestyle='solid', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['LiF_2'][1], data_dict['LiF_1'][1], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['B4C_1'][1], color='red', linestyle='--', label=None)
    plt.plot(thicknesses_in_mm, data_dict['B4C_2'][1], color='red', linestyle='--', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['B4C_2'][1], data_dict['B4C_1'][1], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][1], color='orange', linestyle='-.', label=None)
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], color='orange', linestyle='-.', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)
    plt.xlabel('Thickness (mm)')
    plt.ylabel('$P_{albedo} + P_{transmission}^2P_{backscattering}$ (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('$fom_{back}$')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    fig.savefig('../output/fom_back_at_%.4f_eV.pdf' % energy, bbox_inches='tight')

    # INTERNAL
    fig = plt.figure()
    plt.plot(thicknesses_in_mm, data_dict['Cd'][2], color='green', linestyle='dotted', label='Cd')
    plt.plot(thicknesses_in_mm, data_dict['Gd2O3'][2], color='black', linestyle='solid', label='Gd$_2$O$_3$')

    plt.plot(thicknesses_in_mm, data_dict['LiF_1'][2], color='blue', linestyle='solid', label=None)
    plt.plot(thicknesses_in_mm, data_dict['LiF_2'][2], color='blue', linestyle='solid', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['LiF_2'][2], data_dict['LiF_1'][2], color='blue',
                     label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['B4C_1'][2], color='red', linestyle='--', label=None)
    plt.plot(thicknesses_in_mm, data_dict['B4C_2'][2], color='red', linestyle='--', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['B4C_2'][2], data_dict['B4C_1'][2], color='red',
                     label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_1'][2], color='orange', linestyle='-.', label=None)
    plt.plot(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][2], color='orange', linestyle='-.', label=None)
    plt.fill_between(thicknesses_in_mm, data_dict['epoxyGd2O3_2'][2], data_dict['epoxyGd2O3_1'][2], color='orange',
                     label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)
    plt.xlabel('Thickness (mm)')
    plt.ylabel('$P_{transmission} + P_{albedo}$ (%)')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('$fom_{internal}$')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlim(1e-3, 1e2)
    plt.ylim(1e-4, 1e2)
    #plt.legend(title='Material')
    fig.savefig('../output/fom_internal_at_%.4f_eV.pdf' % energy, bbox_inches='tight')

#get_foms(data_array, 0.082)
#get_foms(data_array, 0.168)
get_foms(data_array, 0.025)
#get_foms(data_array, 1)
print('Wop woop')

# ==============================================================================
#            GET TRANSMISSION VS LAMBDA FOR SINGLE AND DOUBLE DENSITY
# ==============================================================================

def absorption_vs_lambda(data, density_in_n_per_cm3, thickness_in_cm, energies):
    absorption_array = hf.get_absorption_cross(data)
    sigma_abs_in_mm2_vec = []
    for energy in energies:
        sigma_abs_in_mm2_vec.append(hf.get_sigma(absorption_array, energy))
    sigma_abs_in_mm2_vec = np.array(sigma_abs_in_mm2_vec)
    p_absorption = hf.get_probability(sigma_abs_in_mm2_vec, density_in_n_per_cm3, thickness_in_cm)
    return p_absorption

pure_Gd2O3_density_g_per_cm3 = 7.407
pure_Gd2O3_density_n_per_cm3 = float(Gd2O3_data['metadata']['NAtomsPerVolume [1/cm3]'])
layer_thickness_in_cm = 12.5e-4
single_Gd2O3_density_g_per_cm3 = (2.4e-3)/(12.5e-4)
double_Gd2O3_density_g_per_cm3 = (4.8e-3)/(12.5e-4)
fraction_single = single_Gd2O3_density_g_per_cm3/pure_Gd2O3_density_g_per_cm3
fraction_double = double_Gd2O3_density_g_per_cm3/pure_Gd2O3_density_g_per_cm3

# Single density
fig = plt.figure()
energies = np.logspace(-4, -1, 100)
single_density_in_n_per_cm3 = pure_Gd2O3_density_n_per_cm3 * fraction_single
for i in np.arange(1, 11, 1):
    thickness_in_cm = i * layer_thickness_in_cm
    absorption = absorption_vs_lambda(Gd2O3_data, single_density_in_n_per_cm3, thickness_in_cm, energies)
    plt.plot(hf.meV_to_A(1000*energies), absorption*100, label='%d × 12.5 μm layer' %i)

plt.xlabel('Wavelength (Å)')
plt.ylabel('Absorption (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Number of layers', loc=4)
plt.title('Absorption probability - Single density')
plt.xlim(2, 20)
plt.ylim(90, 100)
plt.hlines(99, xmin=2, xmax=20, label=None, color='black', zorder=20)
fig.savefig('../output/absorption_vs_lambda_single_density.pdf', bbox_inches='tight')

# Double density
fig = plt.figure()
energies = np.logspace(-4, -1, 100)
double_density_in_n_per_cm3 = pure_Gd2O3_density_n_per_cm3 * fraction_double
for i in np.arange(1, 11, 1):
    thickness_in_cm = i * layer_thickness_in_cm
    absorption = absorption_vs_lambda(Gd2O3_data, double_density_in_n_per_cm3, thickness_in_cm, energies)
    plt.plot(hf.meV_to_A(1000*energies), absorption*100, label='%d × 12.5 μm layer' %i)

plt.xlabel('Wavelength (Å)')
plt.ylabel('Absorption (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Number of layers', loc=4)
plt.title('Absorption probability - Double density')
plt.xlim(2, 20)
plt.ylim(90, 100)
plt.hlines(99, xmin=2, xmax=20, label=None, color='black', zorder=20)
fig.savefig('../output/absorption_vs_lambda_double_density.pdf', bbox_inches='tight')

print('Single density: %f' % single_Gd2O3_density_g_per_cm3)
print('Double density: %f' % double_Gd2O3_density_g_per_cm3)
print('Ratio 1: %f' % (single_Gd2O3_density_g_per_cm3/pure_Gd2O3_density_g_per_cm3))
print('Ratio 2: %f' % (double_Gd2O3_density_g_per_cm3/pure_Gd2O3_density_g_per_cm3))

energies = np.linspace(2, 200, 100) * 1e-3
thicknesses = np.arange(2, 14, 2) * (1e-6) * (1e2)
fig = plt.figure()
for thickness_in_cm in thicknesses:
    absorption = absorption_vs_lambda(B4C_99_perc_data, B4C_99_perc_density,
                                      thickness_in_cm, energies)
    plt.plot(energies*1000, absorption*100, label='%.2f μm' % (thickness_in_cm*10000))
plt.vlines(82, 0, 100, color='black', label=None)
plt.vlines(167, 0, 100, color='black', label=None)
plt.xlabel('Energy (meV)')
plt.ylabel('Absorption (%)')
plt.xlim(2, 200)
plt.ylim(0, 100)

plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Radial coating', loc=1)
plt.title('Absorption probability')
fig.savefig('../output/absorption_radial_coating.pdf', bbox_inches='tight')

print('DONE!')


print('... and now for T-REX energies')


# Single density
fig = plt.figure()
energies = np.logspace(-2, 0, 100)
single_density_in_n_per_cm3 = pure_Gd2O3_density_n_per_cm3 * fraction_single
for i in np.arange(10, 41, 10):
    thickness_in_cm = i * layer_thickness_in_cm
    absorption = absorption_vs_lambda(Gd2O3_data, single_density_in_n_per_cm3, thickness_in_cm, energies)
    plt.plot(hf.meV_to_A(1000*energies), absorption*100, label='%d × 12.5 μm layer (%d μm)' % (i, int(12.5*i)))

plt.xlabel('Wavelength (Å)')
plt.ylabel('Absorption (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Number of layers', loc=4)
plt.title('Absorption probability - Single density')
plt.xlim(0.5, 1.5)
plt.ylim(30, 100)
plt.hlines(99, xmin=0.5, xmax=6.5, label=None, color='black', zorder=20)
plt.vlines(0.7, ymin=0, ymax=100, label=None, color='black', zorder=20)
fig.savefig('../output/absorption_vs_lambda_single_density_thermal_energies.pdf', bbox_inches='tight')

# Double density
fig = plt.figure()
energies = np.logspace(-2, 0, 100)
double_density_in_n_per_cm3 = pure_Gd2O3_density_n_per_cm3 * fraction_double
for i in np.arange(10, 41, 10):
    thickness_in_cm = i * layer_thickness_in_cm
    absorption = absorption_vs_lambda(Gd2O3_data, double_density_in_n_per_cm3, thickness_in_cm, energies)
    plt.plot(hf.meV_to_A(1000*energies), absorption*100, label='%d × 12.5 μm layer (%d μm)' % (i, int(12.5*i)))

plt.xlabel('Wavelength (Å)')
plt.ylabel('Absorption (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Number of layers', loc=4)
plt.title('Absorption probability - Double density')
plt.xlim(0.5, 1.5)
plt.ylim(30, 100)
plt.hlines(99, xmin=0.5, xmax=6.5, label=None, color='black', zorder=20)
plt.vlines(0.7, ymin=0, ymax=100, label=None, color='black', zorder=20)
fig.savefig('../output/absorption_vs_lambda_double_density_thermal_energies.pdf', bbox_inches='tight')

print('...done')



# ==============================================================================
#            GET TRANSMISSION VS LAMBDA FOR B4C (97% enriched)
# ==============================================================================

fig = plt.figure()
energies_in_eV = np.logspace(-4, -1, 100)
b4c_density_in_n_per_cm3 = float(B4C_97_perc_data['metadata']['NAtomsPerVolume [1/cm3]'])

layer_thickness_in_um = 10
layer_thickness_in_cm = (layer_thickness_in_um) * (1e-6) * (1e2)
for i in np.arange(1, 11, 1):
    thickness_in_cm = i * layer_thickness_in_cm
    absorption = absorption_vs_lambda(B4C_97_perc_data, b4c_density_in_n_per_cm3,
                                      thickness_in_cm, energies_in_eV)
    plt.plot(hf.meV_to_A(1000*energies_in_eV), absorption*100, label='%d μm' % int(layer_thickness_in_um*i))

plt.xlabel('Wavelength (Å)')
plt.ylabel('Absorption (%)')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Thickness', loc=4)
plt.title('Absorption probability - B$_4$C (97 at.% $^{10}$B in B)')
plt.xlim(2, 20)
plt.ylim(90, 100)
plt.hlines(99, xmin=2, xmax=20, label=None, color='black', zorder=20)
fig.savefig('../output/absorption_vs_lambda_b4c_97_percent.pdf', bbox_inches='tight')


# ==============================================================================
#               GET THICKNESS AND ALBEDO FOR A FIXED TRANSMISSION LEVEL
# ==============================================================================

print('Okej då börjar vi!')

transmission_level = 0.01
thicknesses = np.logspace(-6, 5, 1000) * 0.1  # The unit must be in cm
energies = np.logspace(-5, 0, 1000)

data_dict = {'LiF_1': None,
             'LiF_2': None,
             'B4C_1': None,
             'B4C_2': None,
             'Cd': None,
             'Gd2O3': None,
             'epoxyGd2O3_1': None,
             'epoxyGd2O3_2': None,
             'boral': None}
for i, (material, data, density) in enumerate(data_array):
    print('%d/%d' % (i+1, len(data_array)))
    required_thicknesses = []
    albedos = []
    for energy in energies:
        scattering_array = hf.get_scattering_cross(data)
        absorption_array = hf.get_absorption_cross(data)
        thickness, albedo_level, __ = hf.get_transmission(transmission_level, energy,
                                                          scattering_array, absorption_array,
                                                          density, thicknesses)
        required_thicknesses.append(thickness*10) # Make unit to mm
        albedos.append(albedo_level*100)
    data_dict[material] = [np.array(required_thicknesses), np.array(albedos)]

hf.set_thick_labels(12)
fig = plt.figure()
fig.suptitle('Transmission = 1 %')
plt.subplot(1, 2, 1)
plt.plot(energies, data_dict['Cd'][0], color='green', linestyle='dotted', label='Cd')
plt.plot(energies, data_dict['Gd2O3'][0], color='black', linestyle='solid', label='Gd$_2$O$_3$')

plt.plot(energies, data_dict['LiF_1'][0], color='blue', linestyle='solid', label=None)
plt.plot(energies, data_dict['LiF_2'][0], color='blue', linestyle='solid', label=None)
plt.fill_between(energies, data_dict['LiF_2'][0], data_dict['LiF_1'][0], color='blue',
                 label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

plt.plot(energies, data_dict['B4C_1'][0], color='red', linestyle='--', label=None)
plt.plot(energies, data_dict['B4C_2'][0], color='red', linestyle='--', label=None)
plt.fill_between(energies, data_dict['B4C_2'][0], data_dict['B4C_1'][0], color='red',
              label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

plt.plot(energies, data_dict['epoxyGd2O3_1'][0], color='orange', linestyle='-.', label=None)
plt.plot(energies, data_dict['epoxyGd2O3_2'][0], color='orange', linestyle='-.', label=None)
plt.fill_between(energies, data_dict['epoxyGd2O3_2'][0], data_dict['epoxyGd2O3_1'][0], color='orange',
              label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

#plt.plot(energies, data_dict['boral'][0], color='grey', linestyle='-', label='Boral')

plt.xlabel('Energy (eV)')
plt.ylabel('Thickness (mm)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-4, 1)
plt.ylim(1e-3, 1e5)
plt.title('Thickness')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
plt.legend(title='Material', loc=2)
#fig.savefig('../output/thickness_vs_energy_at_%.3f_transmission.pdf' % transmission_level, bbox_inches='tight')

plt.subplot(1, 2, 2)
plt.plot(energies, data_dict['Cd'][1], color='green', linestyle='dotted', label='Cd')
plt.plot(energies, data_dict['Gd2O3'][1], color='black', linestyle='solid', label='Gd$_2$O$_3$')

plt.plot(energies, data_dict['LiF_1'][1], color='blue', linestyle='solid', label=None)
plt.plot(energies, data_dict['LiF_2'][1], color='blue', linestyle='solid', label=None)
plt.fill_between(energies, data_dict['LiF_2'][1], data_dict['LiF_1'][1], color='blue',
                 label='LiF ($^6$Li: 0.08-0.99)', alpha=0.3)

plt.plot(energies, data_dict['B4C_1'][1], color='red', linestyle='--', label=None)
plt.plot(energies, data_dict['B4C_2'][1], color='red', linestyle='--', label=None)
plt.fill_between(energies, data_dict['B4C_2'][1], data_dict['B4C_1'][1], color='red',
              label='B$_4$C ($^{10}$B: 0.20-0.99)', alpha=0.3)

plt.plot(energies, data_dict['epoxyGd2O3_1'][1], color='orange', linestyle='-.', label=None)
plt.plot(energies, data_dict['epoxyGd2O3_2'][1], color='orange', linestyle='-.', label=None)
plt.fill_between(energies, data_dict['epoxyGd2O3_2'][1], data_dict['epoxyGd2O3_1'][1], color='orange',
              label='epoxy-Gd$_2$O$_3$ (0.9-0.1 to 0.1-0.9)', alpha=0.3)

#plt.plot(energies, data_dict['boral'][1], color='grey', linestyle='-', label='Boral')

plt.xlabel('Energy (eV)')
plt.ylabel('Albedo (%)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-4, 1)
plt.ylim(1e-3, 1e2)
plt.title('Albedo')
plt.grid(True, which='major', linestyle='--', zorder=0)
plt.grid(True, which='minor', linestyle='--', zorder=0)
#plt.legend(title='Material', loc=2)

fig.set_figheight(5)
fig.set_figwidth(12)
plt.tight_layout()
fig.savefig('../output/albedo_and_thickness_at_%.4f_transmission.pdf' % transmission_level, bbox_inches='tight')

print('FARDIGT woop woop!')


# ==============================================================================
#       GET THICKNESS AND ALBEDO FOR A FIXED TRANSMISSION LEVEL - NEW
# ==============================================================================

transmission_levels = [0.1, 0.01, 0.001, 0.0001]
thicknesses = np.logspace(-6, 5, 1000) * 0.1  # The unit must be in cm
energies = np.logspace(-5, 0, 1000)
titles = {'B4C_1': 'B$_4$C ($^{10}$B: 0.20)',
          'B4C_2': 'B$_4$C ($^{10}$B: 0.99)',
          'Cd': 'Cd',
          'Gd2O3': 'Gd$_2$O$_3$',
          'LiF_1': 'LiF ($^6$Li: 0.08)',
          'LiF_2': 'LiF ($^6$Li: 0.99)',
          'epoxyGd2O3_1': 'epoxy-Gd$_2$O$_3$ (0.9-0.1)',
          'epoxyGd2O3_2': 'epoxy-Gd$_2$O$_3$ (0.1-0.9)',
          'boral': 'Boral'
          }

data_dict = {'LiF_1': None,
             'LiF_2': None,
             'B4C_1': None,
             'B4C_2': None,
             'Cd': None,
             'Gd2O3': None,
             'epoxyGd2O3_1': None,
             'epoxyGd2O3_2': None,
             'boral': None}

hf.set_thick_labels(14)
colors = ['red', 'orange', 'green', 'blue']
linestyles = ['-', '--', 'dotted', '-.']
percentages = ['10 %', '1 %', '0.1 %', '0.01 %']
for i, (material, data, density) in enumerate(data_array):
    print('%d/%d' % (i+1, len(data_array)))
    fig = plt.figure()
    for j, transmission_level in enumerate(transmission_levels):
        required_thicknesses_lst = []
        for energy in energies:
            scattering_array = hf.get_scattering_cross(data)
            absorption_array = hf.get_absorption_cross(data)
            thickness, __, __ = hf.get_transmission(transmission_level, energy,
                                                    scattering_array, absorption_array,
                                                    density, thicknesses)
            required_thicknesses_lst.append(thickness*10) # Make unit to mm
        required_thicknesses = np.array(required_thicknesses_lst)
        plt.plot(energies, required_thicknesses, color=colors[j],
                 linestyle=linestyles[j],
                 label=percentages[j])

    plt.xlabel('Energy (eV)')
    plt.ylabel('Thickness (mm)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-4, 1)
    plt.ylim(1e-3, 1e2)
    plt.title(titles[material])
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.legend(title='Transmission', loc=0)
    fig.savefig('../output/new_trans_function_of_thickness_and_energy_%s.pdf' % material, bbox_inches='tight')





#print('Thicknesses Cd (transmission):')
#print(data_dict['Cd'][0])
# fig = plt.figure()
# fig.set_figheight(5)
# fig.set_figwidth(10)
# for label, albedo_levels, depths in zip(labels,  albedos_vec, required_thicknesses_vec):
#     plt.subplot(1, 2, 1)
#     plt.plot(energies, depths*10, label=label)
#     plt.subplot(1, 2, 2)
#     plt.plot(energies, albedo_levels, label=label)
#     print(albedo_levels)
# plt.subplot(1, 2, 1)
# plt.xlabel('Energy (eV)')
# plt.ylabel('Thickness (mm)')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(1e-5, 1)
# plt.title('$P_{transmission}$=%.3f' % transmission_level)
# plt.grid(True, which='major', linestyle='--', zorder=0)
# plt.grid(True, which='minor', linestyle='--', zorder=0)
# plt.legend()
# plt.subplot(1, 2, 2)
# plt.xlabel('Energy (eV)')
# plt.ylabel('$P_{albedo}$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(1e-5, 1)
# plt.title('$P_{albedo}$')
# plt.grid(True, which='major', linestyle='--', zorder=0)
# plt.grid(True, which='minor', linestyle='--', zorder=0)
# plt.legend()
# plt.tight_layout()
# fig.savefig('../output/albedo_and_transmission.pdf', bbox_inches='tight')
