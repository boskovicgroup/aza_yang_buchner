#!/usr/bin/env python
# coding: utf-8

# ## Azetidinol Formation Data Formatting
# 
# 0. Imports

# In[1]:


import os
import ipysheet
import ipywidgets as widgets
from kinetics.json_definitions import *
from kinetics.Cary50_spectrum import readSpectraFromFile
from kinetics.util import writeJSONFile


# 1. Parse the spectrum files belonging to the calibration experiments. In this experiment, each standard mixture was prepared with only one of the absorbing components.

# In[2]:


# Specify the directory containing the raw spectrum files.
path_to_calibration_spectra = "azetidinol_formation/calibration_raw/acetophenone_piperidinium_tosylate_08032021.csv"
# Specify the target directory; where JSON files will be written.
path_to_calibration_json = "azetidinol_formation/calibration_json"

# Create an array with the name of each determined component
components = ['piperidine_acetophenone', 'tosic_acid']

# Read spectra
wavelengths, spectra = readSpectraFromFile(path_to_calibration_spectra)

# Create a spreadsheet to input concentration values
labels = list(spectra.keys())
num_rows = len(labels)+1
num_cols = len(components)+1
input_prompt = 'Enter the concentration (uM) of each component in each standard'
concentrations_sheet = ipysheet.sheet(rows=num_rows, columns=num_cols, column_headers=False, row_headers=False)
labels_column = ipysheet.column(0, labels, row_start=1, background_color="silver", read_only=True)
components_row = ipysheet.row(0, components, column_start=1, background_color="silver", read_only=True)
concentration_cells = ipysheet.cell_range([[0]*(num_cols-1)]*(num_rows-1), row_start=1, column_start=1)
submit_button = widgets.Button(
    description='Submit',
    disabled=False,
    button_style='',
    tooltip='Write JSON File',
    icon='check'
)
def on_submit(_):
    for i in range(len(labels)):
        concentrations = list(map(lambda concentration: concentration*1e-6, concentration_cells.value[i]))
        spectrum = list(zip(wavelengths, spectra[labels[i]]))
        spectrum_json = standardSpectrumToJSON(components, concentrations, spectrum)
        path_to_write = os.path.join(path_to_calibration_json, os.path.basename(path_to_calibration_spectra)[:-4] +'_'+ labels[i] + '.json')
        writeJSONFile(spectrum_json, path_to_write)
submit_button.on_click(on_submit)

# #Display spreadsheet
print(input_prompt)
widgets.VBox([concentrations_sheet, submit_button])


# 2. Parse the spectrum files belonging to the kinetic experiments. Each experiment has a directory associated to it that contains the spectrum files obtained at each time point during the experiment.

# In[2]:


# Specify the directory containing the raw spectrum files.
path_to_kinetic_spectra = "azetidinol_formation/kinetic_raw/bandpass_280/concentration_variation/05_75uM_280nm_bp_09012021/05_75uM_280nm_bp_09012021.csv"
# Specify the target directory; where JSON files will be written.
path_to_kinetic_json = "azetidinol_formation/kinetic_json/bandpass_280/concentration_variation/05_75uM_280nm_bp_09012021"

experiment_id = "05_75uM_280nm_bp_09012021"

# Read spectra
wavelengths, spectra = readSpectraFromFile(path_to_kinetic_spectra)

# Create a spreadsheet to input time point values and dilution factors
labels = list(spectra.keys())
num_rows = len(labels)+1
num_cols = 2
input_prompt = 'Enter the dilution factor associated to each experiment'
kinetic_sheet = ipysheet.sheet(rows=num_rows, columns=num_cols, column_headers=False, row_headers=False)
filename_column = ipysheet.column(0, labels, row_start=1, background_color="silver", read_only=True)
first_row = ipysheet.row(0, ['', 'Dilution factor'], column_start=0, background_color="silver", read_only=True)
dilution_factor_cells = ipysheet.cell_range([[1]]*(num_rows-1), row_start=1, column_start=1)
submit_button = widgets.Button(
    description='Submit',
    disabled=False,
    button_style='',
    tooltip='Write JSON File',
    icon='check'
)
def on_submit(_):
    for i in range(len(labels)):
        spectrum = list(zip(wavelengths, spectra[labels[i]]))
        time_point = float(labels[i][2:])
        dilution_factor = dilution_factor_cells.value[i][0]
        spectrum_json = kineticSpectrumToJSON(experiment_id, time_point*60, dilution_factor, spectrum)
        path_to_write = os.path.join(path_to_kinetic_json, os.path.basename(path_to_kinetic_spectra)[:-4] +'_'+ labels[i] + '.json')
        writeJSONFile(spectrum_json, path_to_write)
submit_button.on_click(on_submit)

# #Display spreadsheet
print(input_prompt)
widgets.VBox([kinetic_sheet, submit_button])


# In[ ]:




