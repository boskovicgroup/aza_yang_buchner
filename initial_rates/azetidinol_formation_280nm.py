#!/usr/bin/env python
# coding: utf-8

# ## Azetidinol Formation with 280 nm cut-off filter
# #### Experiment: Rate of reaction from UV-Vis spectra
# 0. Imports

# In[1]:


from kinetics.FSQ import compFTOfSpectrum, compSortedEigenSet, compNumEigenVectors, compSEE
from kinetics.initial_rates import compInitialRate, compLeastSquaresRateLaw
from kinetics.json_definitions import readCalibrationData, readKineticData
import os
import numpy
from collections import defaultdict
from matplotlib import pyplot as plt


# 1. Load the calibration data.

# In[2]:


path_to_calibration_data = "azetidinol_formation/calibration_json"
calibration_data = readCalibrationData(path_to_calibration_data)


# 2. Create a $w$ by $m$ matrix $M$ holding the absorbance measurements at each wavelength $w$ for each standard mixture $m$, and a corresponding $n$ by $m$ matrix $C$ with the concentrations of each one of the $n$ components in each of the $m$ standard mixtures.

# In[3]:


M = numpy.array(calibration_data["absorbances"])
C = numpy.array(calibration_data["concentrations"])
wavelengths = calibration_data["wavelengths"]
component_ids = calibration_data["component_ids"]
# wavelength_offset = 601-567
# M = M[wavelength_offset:,:]


# 3. Compute matrix $F$ by taking the Fourier transform of each of the columns of the absorbance matrix $M$. Select only a consecutive subset of the terms in the Fourier transform by specifying the index of the term at the beginning and the number of terms to include.

# In[4]:


path_length = 1.0 #cm
starting_at_term = 2
number_of_terms = 128

F = [
    compFTOfSpectrum(
        #Fourier terms obtained from real, symmetric spectra
        numpy.hstack((numpy.flip(spectrum), spectrum))[:-1],
    ) 
    for spectrum in M.T
]

F = (numpy.array(F).T)[starting_at_term:(starting_at_term + number_of_terms), :]


# 4. Compute the set of eigenvalues and eigenvectors of the variance-covariance matrix.

# In[5]:


eigenvalues, eigenvectors = compSortedEigenSet(F)


# 5. Compute the number of eigenvectors to use and arrange them in a matrix $V$.

# In[6]:


k = compNumEigenVectors(eigenvalues)
V = eigenvectors[:,:k]
print(k)


# 6. Compute $Z = V^TF$, $P = CZ^T(ZZ^T)^{-1}$ and $M = PV^{T}$.

# In[7]:


Z = numpy.matmul(V.T,F)
P = numpy.matmul(C, numpy.matmul(Z.T,numpy.linalg.inv(numpy.matmul(Z,Z.T))))
M = numpy.matmul(P, V.T)


# 7. Compute the standard error of estimation, $SEE$.

# In[8]:


SEE = compSEE(C,M,F,k)
print(SEE)


# 8. Compute the concentration of each component in the mixture at all time points by reading the data from each kinetic experiment file.

# In[9]:


path_to_kinetic_data = 'azetidinol_formation/kinetic_json/conc_var_280nm'
experiments = {}

# Iterate over the kinetic experiment directories.
for directory_name in [directory for directory in os.listdir(path_to_kinetic_data) if not directory.startswith('.')]:
    path_to_directory = os.path.join(path_to_kinetic_data, directory_name)
    
    # Read kinetic experiment data.
    kinetic_experiment = readKineticData(path_to_directory)
    
    # Read time points as relative to the initial time point.
    time_points = kinetic_experiment["time_points"]
    t0 = min(time_points)
    if t0 > 0:
        time_points = list(numpy.array(time_points)-t0)
    
    # Create a matrix holding the absorbance measurements of the mixture.
    mixture_absorbances = numpy.array(kinetic_experiment["absorbances"])
#     mixture_absorbances = numpy.array(kinetic_experiment["absorbances"])[wavelength_offset:,:]
    
    # Compute the Fourier transform of the mixture absorbance matrix.
    mixture_absorbances_FT = [
        compFTOfSpectrum(
            numpy.hstack((numpy.flip(spectrum), spectrum))[:-1], 
            path_length
        ) for spectrum in mixture_absorbances.T
    ]
    
    # Select the same subset of Fourier terms as the calibration absorbance matrix.
    mixture_absorbances_FT = (numpy.array(mixture_absorbances_FT).T)[
        starting_at_term:(starting_at_term + number_of_terms), :]
    
    # Compute concentrations of each component in the mixture.
    concentrations = numpy.matmul(M, mixture_absorbances_FT).real
    
    # Scale concentrations by the dilution factor.
    dilution_factors = numpy.array(kinetic_experiment["dilution_factors"])
    concentrations = concentrations * numpy.reciprocal(dilution_factors)
    
    # Replace negative concentrations values with 0.
    concentrations[concentrations < 0] = 0
    
    # Record time points and computed component concentrations in the experiments dictionary
    experiment_id = kinetic_experiment["experiment_id"]
    experiments[experiment_id] = {
        "time_points": time_points,
        "concentrations": [list(c) for c in concentrations],
    }


# 9. Compute the initial rates by estimating the initial rate of change in concentration of each monitored component.

# In[10]:


experiments.keys()


# In[13]:


# Add reaction information associated to each component ordered according to its position in the component id list.
component_is_reactant = [True, False]
stoichiometric_coefficients = [1.0, 1.0]
# Set parameters used to compute the initial rate of change in concentration of each component in each experiment.
experiments['05_150uM_08192021']['initial_rate_parameters'] = {
    'piperidine_acetophenone': {
        'time_range': slice(0,len(experiments['05_150uM_08192021']['time_points'])),
        'reaction_percentage': 99
    },
    'tosic_acid': {
        'time_range': slice(0,len(experiments['05_150uM_08192021']['time_points'])),
        'reaction_percentage': 99
    },
}
experiments['01_200uM_08162021']['initial_rate_parameters'] = {
    'piperidine_acetophenone': {
        'time_range': slice(0,len(experiments['01_200uM_08162021']['time_points'])),
        'reaction_percentage': 99
    },
    'tosic_acid': {
        'time_range': slice(0,len(experiments['01_200uM_08162021']['time_points'])),
        'reaction_percentage': 99
    },
}
experiments['02_100uM_08162021']['initial_rate_parameters'] = {
    'piperidine_acetophenone': {
        'time_range': slice(0,len(experiments['02_100uM_08162021']['time_points'])),
        'reaction_percentage': 99
    },
    'tosic_acid': {
        'time_range': slice(0,len(experiments['02_100uM_08162021']['time_points'])),
        'reaction_percentage': 99
    },
}
experiments['04_75uM_08182021']['initial_rate_parameters'] = {
    'piperidine_acetophenone': {
        'time_range': slice(1,len(experiments['04_75uM_08182021']['time_points'])),
        'reaction_percentage': 99
    },
    'tosic_acid': {
        'time_range': slice(0,len(experiments['04_75uM_08182021']['time_points'])),
        'reaction_percentage': 99
    },
}
experiments['03_50uM_08172021']['initial_rate_parameters'] = {
    'piperidine_acetophenone': {
        'time_range': slice(1,len(experiments['03_50uM_08172021']['time_points'])),
        'reaction_percentage': 99
    },
    'tosic_acid': {
        'time_range': slice(0,len(experiments['03_50uM_08172021']['time_points'])),
        'reaction_percentage': 99
    },
}
experiment_ids = list(experiments.keys())
for experiment_id in experiment_ids:
    initial_rates = []
    for component_index, component_id in enumerate(component_ids):
        experiment = experiments[experiment_id]
        time_range = experiment['initial_rate_parameters'][component_id]['time_range']
        reaction_percentage = experiment['initial_rate_parameters'][component_id]['reaction_percentage']
        time_points_subset = experiment['time_points'][time_range]
        concentrations_subset = experiment['concentrations'][component_index][time_range]
        stoichiometric_coefficient = stoichiometric_coefficients[component_index]
        t = time_points_subset[:]
        time_points_subset = list(numpy.array(time_points_subset)-time_points_subset[0])

        initial_rate = compInitialRate(stoichiometric_coefficient, time_points_subset, concentrations_subset, reaction_percentage)
        initial_rates.append(initial_rate)
        
        # Generate initial rate curve.
        is_reactant = component_is_reactant[component_index]
        initial_rate_curve = concentrations_subset[0]  + numpy.multiply(
            (-1 if is_reactant else 1) * stoichiometric_coefficient * initial_rate,
            time_points_subset)

        # Plot kinetic data and the initial rate curve.
        plt.title(f"{experiment_id}") 
        plt.xlabel("t/s") 
        plt.ylabel(f"[{component_id}]/M") 
        plt.plot(t, concentrations_subset, "-bo", linestyle='none', label = f"[{component_id}] at t")
        plt.plot(t, initial_rate_curve, "-r", label = "Initial rate")
        plt.legend(loc = "upper right")
        plt.ylim(min(concentrations_subset), max(concentrations_subset))
        plt.xlim(0, t[-1])
        plt.show()
        print("initial_rate = ", initial_rate, "M/s")
    experiments[experiment_id]['initial_rates'] = initial_rates


# 10. Compute the rate constant and order using linear regression.

# In[12]:


# In this cell we wanted to treat initial rates as in a normal non-photochemical reaction,
# but the orders we kept getting were not interpretable

experiments['01_200uM_08162021']['initial_concentrations'] = [200, 0]
experiments['03_50uM_08172021']['initial_concentrations'] = [50, 0]
experiments['02_100uM_08162021']['initial_concentrations'] = [100, 0]
experiments['04_75uM_08182021']['initial_concentrations'] = [75, 0]
experiments['05_150uM_08192021']['initial_concentrations'] = [150, 0]

initial_rates = [
    experiments[experiment_id]['initial_rates'][0]
    for experiment_id in experiment_ids
]
initial_concentrations = [
    [numpy.divide(experiments[experiment_id]['initial_concentrations'][0], 1e6)]
    for experiment_id in experiment_ids
]

[rate_constant, m_piperidine_acetophenone] = compLeastSquaresRateLaw(initial_concentrations, initial_rates)
fit = [
    numpy.log(rate_constant) + numpy.log(c[0])*m_piperidine_acetophenone
    for c in initial_concentrations
]

print(f"piperidine_acetophenone order: {m_piperidine_acetophenone}")
print(f"reaction rate constant: {rate_constant}")

plt.xlabel("Log (Conc.)") 
plt.ylabel("Log (Initial rate)") 
plt.plot(numpy.log(initial_concentrations), numpy.log(initial_rates), "-bo", linestyle='none')
plt.plot(numpy.log(initial_concentrations), fit, "-r")
plt.show()


# In[ ]:




