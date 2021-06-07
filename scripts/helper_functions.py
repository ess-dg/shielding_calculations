import os
import numpy as np
import Core.Units
import Core.Constants
from decimal import Decimal
import XSectParse.ParseXSectFile as x_parse
import XSectParse.PlotXSectFile as x_plot
import matplotlib.pyplot as plt

def get_transmission(transmission_level, energy, cross_scatt, cross_abs, density, thicknesses):
    # Extract cross-sections
    sigma_scatt = get_sigma(cross_scatt, energy)
    sigma_abs = get_sigma(cross_abs, energy)
    # Get scattering and absorption probabilities
    p_absorption = get_probability(sigma_abs, density, thicknesses)
    p_scattering = get_probability(sigma_scatt, density, thicknesses)
    # Get transmission and albedo probabilities
    p_transmission = (1 - p_absorption) * (1 - 0.5*p_scattering)
    p_albedo = ((1 - p_absorption[1:]) ** 2) * 0.5 * np.diff(p_scattering)
    #print('P_albedo:')
    #print(p_albedo)
    # Get closest value to transmission and corresponding albedo
    idx = get_idx_closest(p_transmission, transmission_level)
    thickness = thicknesses[idx]
    albedo = np.sum(p_albedo[:idx])
    transmission = p_transmission[idx]

    #print('Length: %d' % len(p_transmission))

    #print('Absorption: %f, Scattering: %f' % (p_absorption[idx], p_scattering[idx]))

    return thickness, albedo, transmission

def get_sigma(array, energy):
    return array[1][get_idx_closest(array[0]/Core.Units.eV, energy)]

def get_probability(sigma, density, thickness):
    return 1 - np.exp(-sigma*(1e-2)*density*thickness)

def get_idx_closest(array, value):
    return min(enumerate(array), key=lambda x: abs(x[1]-value))[0]

def get_absorption_cross(data, offset=3):
    # Extract energies
    absorption_energies = np.transpose(np.array(data['procs']['Total']))[0]
    total_energies = np.transpose(np.array(data['procs']['Total']))[0]
    scattering_energies = np.transpose(np.array(data['procs']['hadElastic']))[0][offset:]
    # Extract cross-sections
    total_cross = np.transpose(np.array(data['procs']['Total']))[1]
    scattering_cross = np.transpose(np.array(data['procs']['hadElastic']))[1]
    # Extract absorption
    absorption_cross = total_cross - scattering_cross[offset:]
    absorption = np.array([absorption_energies, absorption_cross])
    # Just cross-check that I'm not doing anything stupid...
    if (sum(total_energies - scattering_energies)) > 0:
        print('Might be something wrong with array operations...')
    return absorption

def get_scattering_cross(data):
    return np.transpose(np.array(data['procs']['hadElastic']))

def meV_to_A(energy):
    return np.sqrt(81.81/energy)

def A_to_meV(wavelength):
    return (81.81 /(wavelength ** 2))

def get_number_density(density_g_per_cm3, molar_mass_g_per_mol, atoms_per_molecule):
    number_density_per_cm3 = (density_g_per_cm3/molar_mass_g_per_mol) * Core.Constants.Avogadro * atoms_per_molecule
    return number_density_per_cm3

def plot_cross_area(array_1, array_2, label, color, linestyle, ax=None):
    if ax is None:
        plt.plot((array_1[0]/Core.Units.eV), array_1[1]/Core.Units.barn,
                 color=color, linestyle=linestyle, label=None)
        plt.plot((array_2[0]/Core.Units.eV), array_2[1]/Core.Units.barn,
                 color=color, linestyle=linestyle, label=None)
        plt.fill_between(array_1[0]/Core.Units.eV, array_2[1]/Core.Units.barn,
                         array_1[1]/Core.Units.barn, color=color, label=label,
                         alpha=0.3)
    else:
        ax.plot((array_1[0]/Core.Units.eV), array_1[1]/Core.Units.barn,
                color=color, linestyle=linestyle, label=None)
        ax.plot((array_2[0]/Core.Units.eV), array_2[1]/Core.Units.barn,
                color=color, linestyle=linestyle, label=None)
        ax.fill_between(array_1[0]/Core.Units.eV, array_2[1]/Core.Units.barn,
                        array_1[1]/Core.Units.barn, color=color, label=label,
                        alpha=0.3)

def plot_prob_area(array_1, array_2, label, color, linestyle, d_1, d_2, thickness):
    prob_1 = 1 - np.exp(-(array_1[1]*(1e-2)*d_1*thickness))
    prob_2 = 1 - np.exp(-(array_2[1]*(1e-2)*d_2*thickness))
    plt.plot((array_1[0]/Core.Units.eV), prob_1,
             color=color, linestyle=linestyle, label=None)
    plt.plot((array_2[0]/Core.Units.eV), prob_2,
             color=color, linestyle=linestyle, label=None)
    plt.fill_between(array_1[0]/Core.Units.eV, prob_2, prob_1, color=color,
                     label=label, alpha=0.3)

def plot_prob_area_thickness(cross_1, cross_2, label, color, linestyle, d_1, d_2, thicknesses):
    prob_1 = 1 - np.exp(-(cross_1*(1e-2)*d_1*thicknesses))
    prob_2 = 1 - np.exp(-(cross_2*(1e-2)*d_2*thicknesses))
    plt.plot(thicknesses*10, prob_1,
             color=color, linestyle=linestyle, label=None)
    plt.plot(thicknesses*10, prob_2,
             color=color, linestyle=linestyle, label=None)
    plt.fill_between(thicknesses*10, prob_2, prob_1, color=color,
                     label=label, alpha=0.3)

def plot_thick_vs_energy_area(array_1, array_2, label, color, linestyle,
                              density_1, density_2, energy_start, energy_stop,
                              capture_probability):
    # Extract cross-sections
    idx_start = min(enumerate(array_1[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_start))[0]
    idx_stop = min(enumerate(array_1[0]/Core.Units.eV), key=lambda x: abs(x[1]-energy_stop))[0]
    cross_sections_1 = array_1[1][idx_start:idx_stop] * (1e-2)
    cross_sections_2 = array_2[1][idx_start:idx_stop] * (1e-2)
    energies = array_1[0][idx_start:idx_stop]/Core.Units.eV
    # Extract thicknesses
    thicknesses_1 = (- np.log(1 - capture_probability) / (cross_sections_1 * density_1)) * 10
    thicknesses_2 = (- np.log(1 - capture_probability) / (cross_sections_2 * density_2)) * 10
    # Plot
    plt.plot(energies, thicknesses_1, color=color, linestyle=linestyle, label=None)
    plt.plot(energies, thicknesses_2, color=color, linestyle=linestyle, label=None)
    plt.fill_between(energies, thicknesses_2, thicknesses_1, color=color,
                     label=label, alpha=0.3)

def set_thick_labels(thickness):
    # Customize matplotlib font sizes
    plt.rc('font', size=thickness)          # controls default text sizes
    plt.rc('axes', titlesize=thickness)     # fontsize of the axes title
    plt.rc('axes', labelsize=thickness)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=thickness)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=thickness)    # fontsize of the tick labels
    plt.rc('legend', fontsize=thickness)    # legend fontsize
    plt.rc('figure', titlesize=thickness)  # fontsize of the figure title
