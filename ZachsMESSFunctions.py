
# Functions I made to improve MESS
def saveData(self):
        
    print('Saving data...')

# Define the path to save the file
    save_path = self.SavePath_edit.text()  # Get the save path from the GUI
    if not save_path:
        print("Save path is not set. Please choose a save path.")
        return

    event_date = self.event_date if self.event_date is not None else "Wavelength"
    event_time = self.event_time if self.event_time is not None else "Intensity"


    file_name = os.path.join(save_path, f"{event_date}_{event_time}.txt")  # Name of the file

    try:
        # Open the file and write the data
        with open(file_name, 'w') as f:
            f.write("Wavelength (nm)\tIntensity\n")  # Add a header
            for w in range(len(self.spectrumX)):
                f.write(f"{self.spectrumX[w]:.6f}\t{self.spectrumY[w]:.6f}\n")
        print(f"Data saved to {file_name}")
    except Exception as e:
        print(f"Failed to save data: {e}")

    print("Data Saved!")



def saveElementTemperatures(self):
        """Save Tlo and Thi for each element to a dictionary."""
        self.element_temperatures = {}

        for element in self.spectral.elemdata.els:
            element_name = getattr(element, 'elementcode', 'Unknown')
            Tlo = getattr(element, 'Tlo', None)
            Thi = getattr(element, 'Thi', None)
            self.element_temperatures[element_name] = (Tlo, Thi)

        print("Element temperatures saved:", self.element_temperatures)
    
def getFittingElements(self):  # Retrieve and view the fitted and locked elements
        print("Debug: Retrieving fitted and locked elements")
        fitting_elements = []
        max_iterations = 1000  # Limit the number of iterations for debugging / KEEP HIGH, will not add important fitting elements if too low
        iteration_count = 0

        try:
            for element in self.spectral.elemdata.els: # Attempting to cap the for loop, could be what is crashing MESS
                iteration_count += 1
                if iteration_count > max_iterations:
                    print("Debug: Max iterations reached, breaking loop.")
                    break

                if not hasattr(element, 'user_fitflag'):
                    print("Warning: Element missing 'user_fitflag' attribute.")
                    continue

                # Include both fitting (1) and locked (2) elements
                if element.user_fitflag in [1, 2]:
                    fitting_elements.append({
                        'name': getattr(element, 'elementcode', 'Unknown'),  # Element name or code
                        'N_warm': getattr(element, 'N_warm', None),  # Warm plasma density
                        'N_hot': getattr(element, 'N_hot', None),    # Hot plasma density
                        'speclo': getattr(element, 'speclo', None),  # Low-temperature spectrum
                        'spechi': getattr(element, 'spechi', None),  # High-temperature spectrum
                        'state': 'locked' if element.user_fitflag == 2 else 'fitting'  # Add state info
                    })
        except Exception as e:
            print(f"Error processing elements in getFittingElements: {e}")
        return fitting_elements


def displayFittingElements(self):  # Puts fitted elements to the output
    try:
        fitting_elements = self.getFittingElements()
        if not fitting_elements:
            print("No fitting elements found.")
            return

        plt.figure()
        for elem in fitting_elements:
            if elem['speclo'] is not None and elem['spechi'] is not None:
                # Convert speclo and spechi from pointers to NumPy arrays
                speclo_array = np.ctypeslib.as_array(elem['speclo'], shape=(self.spectral.spcalib.nwavelengths,))
                spechi_array = np.ctypeslib.as_array(elem['spechi'], shape=(self.spectral.spcalib.nwavelengths,))
                
                # Plot the data
                plt.plot(speclo_array, spechi_array, label=elem['name'])

        if fitting_elements:
            plt.legend()
        else:
            print("No elements to display in the legend.")

            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Intensity")
            plt.title("Fitting Elements")
            plt.draw()  # Use draw instead of show to avoid event loop conflicts
    except Exception as e:
        print(f"Error displaying fitting elements: {e}")

# Attempting to save the fitted elements
def saveFittedElements(self):  # Save fitted and locked elements along with temperature parameters
    try:
        # Retrieve the fitted elements
        fitting_elements = self.getFittingElements()
        if not fitting_elements:
            print("No fitting elements found.")
            return

        # Define the path to save the file
        save_path = self.SavePath_edit.text()  # Get the save path from the GUI
        if not save_path:
            print("Save path is not set. Please choose a save path.")
            return

        file_name = os.path.join(save_path, "fitted_elements.txt")

        max_spec_points = 10  # Limit the number of spectral points to save

        # Write the fitted elements and temperature parameters to a file
        with open(file_name, 'w') as f:
            f.write("Fitted Elements Data\n")
            f.write("====================\n")
            # Saving temperature parameters at end of fit
            f.write(f"Warm Temperature (Tlo): {self.spectral.elemdata.Tlo} K\n")
            f.write(f"Hot Temperature (Thi): {self.spectral.elemdata.Thi} K\n")
            f.write("\n")

            for elem in fitting_elements:
                f.write(f"Element: {elem['name']}\n")
                f.write(f"  State: {elem['state']}\n")
                f.write(f"  N_warm: {elem['N_warm']}\n")
                f.write(f"  N_hot: {elem['N_hot']}\n")

                # Converting speclo and spechi from pointers to NumPy arrays
                if elem['speclo'] is not None and elem['spechi'] is not None:
                    try:
                        speclo_array = np.ctypeslib.as_array(elem['speclo'], shape=(self.spectral.spcalib.nwavelengths,))
                        spechi_array = np.ctypeslib.as_array(elem['spechi'], shape=(self.spectral.spcalib.nwavelengths,))
                        if len(speclo_array) > max_spec_points:
                            f.write("  (Only first 10 points shown)\n")
                            f.write("  speclo: " + ", ".join(f"{val:.6f}" for val in speclo_array[:max_spec_points]) + "\n")
                            f.write("  spechi: " + ", ".join(f"{val:.6f}" for val in spechi_array[:max_spec_points]) + "\n")
                    except Exception as e:
                        f.write(f"  Error converting speclo/spechi: {e}\n")
                else:
                    f.write("  speclo: None\n")
                    f.write("  spechi: None\n")

                f.write("\n")

        print(f"Fitted and locked elements saved to {file_name}")
    except Exception as e:
        print(f"Error saving fitted and locked elements: {e}")

# Calculate the average intensity of Fe lines in specified ranges
def calculateAverageFeIntensity(self, fe_ranges=[(430, 440), (490, 515), (520, 550)]):
    """
    Calculate the average intensity of Fe lines within a given wavelength range.

    Arguments:
        fe_ranges: [list of tuples] Wavelength ranges for Fe lines (e.g., [(500, 520), (430, 440)]).

    Returns:
        fe_averages: [dict] Average intensities for each range.
    """
    fe_ranges = [(430, 440), (490, 515), (520, 550)]

    fe_averages = {}  # Store average intensities for each range

    # Validate fe_ranges
    if not isinstance(fe_ranges, list) or not all(isinstance(r, tuple) and len(r) == 2 for r in fe_ranges):
        print(f"Invalid fe_ranges argument: {fe_ranges}. Expected a list of tuples.")
        return None

    try:
        # Extract the spectrum data
        wavelengths = np.array(self.spectrumX)
        intensities = np.array(self.spectrumY_resp)

        # Filter the spectrum for the Fe line region
        for fe_range in fe_ranges:
            print(f"Calculating average intensity for Fe lines in range: {fe_range}")
            fe_mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
            fe_intensities = intensities[fe_mask]

            # Calculate the average intensity
            average_intensity = np.mean(fe_intensities) if len(fe_intensities) > 0 else 0.0
            if float(average_intensity) == 0.0:
                print(f"No intensities found in range {fe_range}. Skipping this range.")
                continue  # Skip if no intensities in this range

            # Store the average intensity in the dictionary
            else:
                fe_averages[fe_range] = average_intensity

            print(f"Average intensity of Fe lines in range {fe_range}: {average_intensity:.6f}")

        return fe_averages
    except Exception as e:
        print(f"Error calculating average Fe intensity: {e}")
        return None



def subtractContinuum(self):
    """
    Perform linear continuum subtraction around Mg (518nm) and Na (589nm).
    """
    try:
        # Define wavelength ranges for continuum subtraction
        mg_continuum_left = (510, 515)  # Left continuum range for Mg
        mg_continuum_right = (520, 525)  # Right continuum range for Mg
        na_continuum_left = (580, 585)  # Left continuum range for Na
        na_continuum_right = (593, 598)  # Right continuum range for Na

        # Calculate average intensity of Fe lines
        fe_ranges=[(430, 440), (490, 515), (520, 550)]
        fe_averages = self.calculateAverageFeIntensity(fe_ranges)

        # Access values from the fe_averages dictionary
        for fe_range, average_intensity in fe_averages.items():
            print(f"Fe range {fe_range}: Average intensity = {average_intensity:.6f}")

        # Calculate the average of Fe intensities
        fe_average_value = np.mean(list(fe_averages.values())) if fe_averages else 0.0
        print(f"Average Fe intensity: {fe_average_value:.6f}")

        # Extract the spectrum data
        wavelengths = np.array(self.spectrumX)
        intensities = np.array(self.spectrumY_resp) 

        # Fit linear continuum for Mg
        mg_left_mask = (wavelengths >= mg_continuum_left[0]) & (wavelengths <= mg_continuum_left[1])
        mg_right_mask = (wavelengths >= mg_continuum_right[0]) & (wavelengths <= mg_continuum_right[1])
        mg_continuum_wavelengths = np.concatenate([wavelengths[mg_left_mask], wavelengths[mg_right_mask]])
        mg_continuum_intensities = np.concatenate([intensities[mg_left_mask], intensities[mg_right_mask]])
        mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
        mg_baseline = np.polyval(mg_fit, wavelengths)

        # Fit linear continuum for Na
        na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
        na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
        na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
        na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
        na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
        na_baseline = np.polyval(na_fit, wavelengths)

        # Subtract the continuum
        corrected_intensities = intensities.copy()
        corrected_intensities[(wavelengths >= 510) & (wavelengths <= 525)] -= mg_baseline[(wavelengths >= 510) & (wavelengths <= 525)]
        corrected_intensities[(wavelengths >= 580) & (wavelengths <= 598)] -= na_baseline[(wavelengths >= 580) & (wavelengths <= 598)]

        # Store the intensities at 518 nm (Mg) and 589 nm (Na) in a dictionary
        intensity_dict = {}
        mg_index = np.argmin(np.abs(wavelengths - 518))  # Find the index closest to 518 nm
        na_index = np.argmin(np.abs(wavelengths - 589))  # Find the index closest to 589 nm
        intensity_dict['Mg_518nm'] = round(corrected_intensities[mg_index], 6)  # Round to 6 decimal places
        intensity_dict['Na_589nm'] = round(corrected_intensities[na_index], 6)

        print(f"Intensity at 518 nm (Mg): {intensity_dict['Mg_518nm']:.6f}")
        print(f"Intensity at 589 nm (Na): {intensity_dict['Na_589nm']:.6f}")

        # Update the spectrum
        self.spectrumY_resp = corrected_intensities

        # Plot the corrected spectrum
        self.clearSpec()
        self.Plot.plot(wavelengths, corrected_intensities, pen=pg.mkPen(color=(0, 0, 255), width=2))
        self.Plot.setLabel('left', 'Intensity')
        self.Plot.setLabel('bottom', 'Wavelength (nm)')
        self.Plot.setXRange(np.min(wavelengths), np.max(wavelengths))
        self.Plot.setYRange(0, np.max(corrected_intensities) * 1.2)

        print("Continuum subtraction completed successfully.")
        print("Intensity dictionary:", intensity_dict)

        intensity_sum = intensity_dict['Mg_518nm'] + intensity_dict['Na_589nm'] + fe_average_value
        print(f"Sum of intensities (Mg + Na + Fe): {intensity_sum:.6f}")

        # With these corrected intensities, we may find normalized values of the elements, which we need for the ternary plots

        norm_intensities = {}

        norm_intensities['Mg'] = round(intensity_dict['Mg_518nm'] / intensity_sum, 6) if intensity_sum != 0 else 0.0
        norm_intensities['Na'] = round(intensity_dict['Na_589nm'] / intensity_sum, 6) if intensity_sum != 0 else 0.0
        norm_intensities['Fe'] = round(fe_average_value / intensity_sum, 6) if intensity_sum != 0 else 0.0
        print(f"Normalized intensities: {norm_intensities}")
        # Store the normalized intensities in the spectral object
        self.spectral.normalized_intensities = norm_intensities

    except Exception as e:
        print(f"Error during continuum subtraction: {e}")



# Put these into the Ui Class

 # Save fitted elements
        self.SaveFittedElements_button.clicked.connect(self.saveFittedElements)
        #MJM
        self.ChooseSavePath_button.clicked.connect(self.chooseSavePath)

        # Fe Average Button
        self.AverageFe_button.clicked.connect(self.calculateAverageFeIntensity)

        # Continuum subtraction
        self.ContinuumSubtract_button.clicked.connect(self.subtractContinuum)

# initialize the GuralSpectral object
        self.spectral = spectral_library.GuralSpectral(10000, 4500, None, None, None, None)      


    # Need to add buttons for the above Ui functions into the UI as well  
