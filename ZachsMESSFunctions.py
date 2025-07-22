
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
    def calculateAverageFeIntensity(self, fe_ranges=[(430, 440), (485, 515), (520, 550)], top_n=15, exclusion_nm=1.0):
        """
        Calculate the average of the highest N intensity points across all Fe ranges,
        excluding any point within ±exclusion_nm of a previously selected peak.
        Returns both the per-range averages and the top-N average.
        """
        fe_ranges = [(420, 455), (485, 515), (520, 550)]
        fe_averages = {}  # Store average intensities for each range
        all_fe_points = []  # Store (wavelength, intensity) tuples

        # Validate fe_ranges
        if not isinstance(fe_ranges, list) or not all(isinstance(r, tuple) and len(r) == 2 for r in fe_ranges):
            print(f"Invalid fe_ranges argument: {fe_ranges}. Expected a list of tuples.")
            return None, [], 0.0

        try:
            wavelengths = np.array(self.spectrumX)
            intensities = np.array(self.spectrumY_resp)

            for fe_range in fe_ranges:
                fe_mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
                fe_waves = wavelengths[fe_mask]
                fe_intensities = intensities[fe_mask]
                all_fe_points.extend(list(zip(fe_waves, fe_intensities)))
                average_intensity = np.mean(fe_intensities) if len(fe_intensities) > 0 else 0.0
                if float(average_intensity) == 0.0:
                    print(f"No intensities found in range {fe_range}. Skipping this range.")
                    continue
                fe_averages[fe_range] = average_intensity
                print(f"Average intensity of Fe lines in range {fe_range}: {average_intensity:.6f}")

            # Sort all points by intensity (descending)
            # all_fe_points = sorted(all_fe_points, key=lambda x: x[1], reverse=True)

            all_fe_points_desc = sorted(all_fe_points, key=lambda x: x[1], reverse=True)
            all_fe_points_asc = sorted(all_fe_points, key=lambda x: x[1])


            # Select the top N points, excluding those within exclusion_nm of any already selected peak
            # selected_peaks = []
            # selected_wavelengths = []

            # for wave, inten in all_fe_points:
            #     # Exclude if within exclusion_nm of any already selected peak
            #     if any(abs(wave - sel_wave) < exclusion_nm for sel_wave in selected_wavelengths):
            #         continue
            #     selected_peaks.append(inten)
            #     selected_wavelengths.append(wave)
            #     if len(selected_peaks) == top_n:
            #         break

            # Select top N (highest) with exclusion
            selected_top_peaks = []
            selected_top_wavelengths = []
            for wave, inten in all_fe_points_desc:
                if any(abs(wave - sel_wave) < exclusion_nm for sel_wave in selected_top_wavelengths):
                    continue
                selected_top_peaks.append(inten)
                selected_top_wavelengths.append(wave)
                if len(selected_top_peaks) == top_n:
                    break

            # Select bottom N (lowest) with exclusion
            selected_bottom_peaks = []
            selected_bottom_wavelengths = []
            for wave, inten in all_fe_points_asc:
                if any(abs(wave - sel_wave) < exclusion_nm for sel_wave in selected_bottom_wavelengths):
                    continue
                selected_bottom_peaks.append(inten)
                selected_bottom_wavelengths.append(wave)
                if len(selected_bottom_peaks) == top_n:
                    break

            top_points = selected_top_peaks
            bottom_points = selected_bottom_peaks
            top_n_average = np.mean(top_points) if top_points else 0.0
            bottom_n_average = np.mean(bottom_points) if bottom_points else 0.0

            print(f"Top {top_n} Fe intensities (with exclusion): {top_points}")
            print(f"Average of top {top_n} Fe intensities: {top_n_average:.6f}")
            print(f"Bottom {top_n} Fe intensities (with exclusion): {bottom_points}")
            print(f"Average of bottom {top_n} Fe intensities: {bottom_n_average:.6f}")


            return all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average
        
        except Exception as e:
            print(f"Error calculating average Fe intensity: {e}")
            return None, [], 0.0



#Keep this function for continuum subtraction
    def subtractContinuum(self):
        """
        Perform linear continuum subtraction around Mg (518nm) and Na (589nm).
        """
        try:
            # Define wavelength ranges for continuum subtraction
            mg_continuum_left = (485, 515)  # Left continuum range for Mg
            mg_continuum_right = (520, 550)  # Right continuum range for Mg
            na_continuum_left = (570, 585)  # Left continuum range for Na
            na_continuum_right = (593, 608)  # Right continuum range for Na

            # Accessing user-specified regions to subtract from
            # mg_left = self.MgLeft_spin.value()
            # mg_right = self.MgRight_spin.value()
            # na_left = self.NaLeft_spin.value()
            # na_right = self.NaRight_spin.value()

            # Calculate average intensity of Fe lines
            fe_ranges=[(420, 455), (485, 515), (520, 550)]
            all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
            
            fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0
            # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0
            # fe_average_value = (top_n_average + bottom_n_average) / 2 if top_n_average and bottom_n_average else 0.0
            print(f"Top N Fe average intensity: {top_n_average:.6f}")
            print(f"Bottom N Fe average intensity: {bottom_n_average:.6f}")
            print(f"Overall Fe average intensity: {fe_average_value:.6f}\n")

            # Get the highest Fe intensity in all Fe ranges
            max_fe_intensity = max([inten for wave, inten in all_fe_points]) if all_fe_points else 0.0
            print(f"Maximum Fe intensity in ranges: {max_fe_intensity:.6f}")

            # Access values from the fe_averages dictionary
            for fe_range, average_intensity in fe_averages.items():
                print(f"Fe range {fe_range}: Average intensity = {average_intensity:.6f}")

            # Calculate the average of Fe intensities
            #fe_average_value = np.mean(list(fe_averages.values())) if fe_averages else 0.0
            print(f"Average Fe intensity: {fe_average_value:.6f}")

            # Extract the spectrum data
            wavelengths = np.array(self.spectrumX)
            intensities = np.array(self.spectrumY_resp) 

            # Integrated Iron Lines
            # fe_lines = [427.3, 430.8, 432.6, 438.4, 440.5, 492.0, 495.7, 504.7, 526.9, 532.8, 537.1, 540.4, 543.1, 544.9]
            # fe_int_before = sum(self.integrate_line(wavelengths, intensities, line, width=1.0) for line in fe_lines)
            # print(f"Integrated Fe intensity before subtraction: {fe_int_before:.6f}")

            # ///////////////////////////////////
            fe_ranges = [(420, 455), (485, 515), (520, 550)]
            fe_int_before = 0.0
            for fe_range in fe_ranges:
                mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
                fe_int_before += np.trapz(intensities[mask], wavelengths[mask])
            print(f"TEST:   Integrated Fe intensity over ranges before subtraction: {fe_int_before:.6f}")

            # METHOD 1: For Comparing Regions

            fe_ranges = [(420, 455), (485, 515), (520, 550)]
            fe_ints = []
            for fe_range in fe_ranges:
                mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
                fe_int = np.trapz(intensities[mask], wavelengths[mask])
                fe_ints.append(fe_int)
            average_fe_int = np.mean(fe_ints)
            print(f"TEST1:   Average integrated Fe intensity per range: {average_fe_int:.6f}")
            print(f"TEST1:   Average integrated Fe intensity per line:  {average_fe_int/15:.6f}\n")

            # METHOD 2: For Normalization
            fe_areas = []
            for fe_range in fe_ranges:
                mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
                area = np.trapz(intensities[mask], wavelengths[mask])
                width = fe_range[1] - fe_range[0]
                fe_areas.append(area / width)
            average_fe_per_nm = np.mean(fe_areas)
            print(f"TEST2:   Average Fe intensity per nm: {average_fe_per_nm:.6f}\n")

            # METHOD 3: For pointwise normalization
            fe_points = []
            for fe_range in fe_ranges:
                mask = (wavelengths >= fe_range[0]) & (wavelengths <= fe_range[1])
                fe_points.extend(intensities[mask])
            average_fe_point = np.mean(fe_points)
            print(f"TEST3:   Average Fe intensity (all points in ranges): {average_fe_point:.6f}\n")

            # ///////////////////////////////////

            # Fit linear continuum for Mg
            mg_left_mask = (wavelengths >= mg_continuum_left[0]) & (wavelengths <= mg_continuum_left[1])
            mg_right_mask = (wavelengths >= mg_continuum_right[0]) & (wavelengths <= mg_continuum_right[1])
            mg_continuum_wavelengths = np.concatenate([wavelengths[mg_left_mask], wavelengths[mg_right_mask]])
            mg_continuum_intensities = np.concatenate([intensities[mg_left_mask], intensities[mg_right_mask]])
            # # For trendline
            # mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
            # mg_baseline = np.polyval(mg_fit, wavelengths)

            # For horizontal line
            mg_level = np.median(mg_continuum_intensities)
            mg_baseline = np.full_like(wavelengths, mg_level)
            

            # Fit linear continuum for Na
            na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
            na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
            na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
            na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
            # # For trendline
            # na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
            # na_baseline = np.polyval(na_fit, wavelengths)

            # For horizontal line
            na_level = np.median(na_continuum_intensities)
            na_baseline = np.full_like(wavelengths, na_level)

            # Subtract the continuum
            corrected_intensities = intensities.copy()
            corrected_intensities[(wavelengths >= 485) & (wavelengths <= 550)] -= mg_baseline[(wavelengths >= 485) & (wavelengths <= 550)]
            corrected_intensities[(wavelengths >= 570) & (wavelengths <= 608)] -= na_baseline[(wavelengths >= 570) & (wavelengths <= 608)]

            # ///////////////////////////////////
            mg_range = (517.0, 519.0)
            na_range = (588.0, 590.0)

            mg_mask = (wavelengths >= mg_range[0]) & (wavelengths <= mg_range[1])
            na_mask = (wavelengths >= na_range[0]) & (wavelengths <= na_range[1])

            print("Mg mask points:", np.sum(mg_mask), "Na mask points:", np.sum(na_mask))
            print("Mg wavelengths:", wavelengths[mg_mask])
            print("Na wavelengths:", wavelengths[na_mask])
            print("Corrected Mg values:", corrected_intensities[mg_mask])
            print("Corrected Na values:", corrected_intensities[na_mask])
            print("Wavelength range:", np.min(wavelengths), np.max(wavelengths))

            mg_vals = np.clip(corrected_intensities[mg_mask], 0, None)
            na_vals = np.clip(corrected_intensities[na_mask], 0, None)

            mg_int = np.trapz(mg_vals, wavelengths[mg_mask])
            na_int = np.trapz(na_vals, wavelengths[na_mask])

            # mg_int = self.integrate_line(wavelengths, corrected_intensities, 518.0, width=1.0)
            # na_int = self.integrate_line(wavelengths, corrected_intensities, 589.0, width=1.0)
            

            print(f"TEST:   Integrated Mg intensity over ranges after subtraction: {mg_int:.6f}")
            print(f"TEST:   Integrated Na intensity over ranges after subtraction: {na_int:.6f}")
            # ///////////////////////////////////

            print("Corrected Mg region:", corrected_intensities[(wavelengths >= 517.5) & (wavelengths <= 518.5)])
            print("Corrected Na region:", corrected_intensities[(wavelengths >= 588.5) & (wavelengths <= 589.5)])

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
            intensity_sum_2 = intensity_dict['Mg_518nm'] + intensity_dict['Na_589nm'] + max_fe_intensity

            print(f"Sum of intensities (Mg + Na + Fe): {intensity_sum_2:.6f}")

            # With these corrected intensities, we may find normalized values of the elements, which we need for the ternary plots

            norm_intensities = {}

            norm_intensities['Fe Max'] = round(max_fe_intensity / intensity_sum_2, 6) if intensity_sum != 0 else 0.0

            # norm_intensities['Fe'] = round(fe_average_value / intensity_sum, 6) if intensity_sum != 0 else 0.0
            norm_intensities['Mg'] = round(intensity_dict['Mg_518nm'] / intensity_sum_2, 6) if intensity_sum != 0 else 0.0
            norm_intensities['Na'] = round(intensity_dict['Na_589nm'] / intensity_sum_2, 6) if intensity_sum != 0 else 0.0
            
            #print(f"Normalized intensities: {norm_intensities}\n")
            print(F"Nromalized intensities (using max fe): {norm_intensities}")
            # Store the normalized intensities in the spectral object
            self.spectral.normalized_intensities = norm_intensities

            # At the end, instead of just printing, return the results:
            return intensity_dict, fe_average_value, norm_intensities

        except Exception as e:
            print(f"Error during continuum subtraction: {e}")
            return None, 0.0, {}

    
    def showContinuumModel(self):
        """
        Displaying the linear models for continuum subtraction around Mg (518nm) and Na (589nm).
        """
        try:
            # Define wavelength ranges for continuum subtraction
            mg_continuum_left = (485, 515)  # Left continuum range for Mg
            mg_continuum_right = (520, 550)  # Right continuum range for Mg
            na_continuum_left = (570, 585)  # Left continuum range for Na
            na_continuum_right = (593, 608)  # Right continuum range for Na

            # Calculate average intensity of Fe lines
            fe_ranges=[(420, 455), (485, 515), (520, 550)]
            all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
            # fe_average_value = (top_n_average + bottom_n_average) / 2 
            # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0
            fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0

            all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
            # # Use bottom_n_average for the Fe trendline
            # fe_trend_y = bottom_n_average
            # fe_trend_x = np.linspace(fe_ranges[0][0], fe_ranges[-1][1], 100)
            # self.Plot.plot(fe_trend_x, [fe_trend_y]*len(fe_trend_x), pen=pg.mkPen(color=(30,144,255), width=2, style=QtCore.Qt.DashLine), name="Fe Trend (Lowest Points)")

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
            # # For trendline
            # mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
            # mg_baseline = np.polyval(mg_fit, wavelengths)

            # Mg horizontal background line
            mg_level = np.median(mg_continuum_intensities)
            mg_baseline = np.full_like(wavelengths, mg_level)
            mg_region_mask = (wavelengths >= 485) & (wavelengths <= 550)
            self.Plot.plot(
                wavelengths[mg_region_mask],
                np.full(np.sum(mg_region_mask), mg_level),
                pen=pg.mkPen(color=(255, 0, 0), width=2, style=QtCore.Qt.DotLine),
                name="Mg Background (mean)"
            )

            # Fit linear continuum for Na
            na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
            na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
            na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
            na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
            # # For trendline
            # na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
            # na_baseline = np.polyval(na_fit, wavelengths)

            # Na horizontal background line
            na_level = np.median(na_continuum_intensities)
            na_baseline = np.full_like(wavelengths, na_level)
            na_region_mask = (wavelengths >= 570) & (wavelengths <= 608)
            self.Plot.plot(
                wavelengths[na_region_mask],
                np.full(np.sum(na_region_mask), na_level),
                pen=pg.mkPen(color=(255, 165, 0), width=2, style=QtCore.Qt.DotLine),
                name="Na Background (mean)"
            )

            # Plot the original spectrum
            self.clearSpec()
            self.Plot.plot(wavelengths, intensities, pen=pg.mkPen(color=(0, 0, 255), width=2), name="Original Spectrum")

            # Plot the Mg continuum model over the Mg region
            mg_region_mask = (wavelengths >= 485) & (wavelengths <= 550)
            self.Plot.plot(wavelengths[mg_region_mask], mg_baseline[mg_region_mask], pen=pg.mkPen(color=(255, 0, 0), width=2, style=QtCore.Qt.DashLine), name="Mg Continuum")

            # Plot the Na continuum model over the Na region
            na_region_mask = (wavelengths >= 570) & (wavelengths <= 608)
            self.Plot.plot(wavelengths[na_region_mask], na_baseline[na_region_mask], pen=pg.mkPen(color=(255, 165, 0), width=2, style=QtCore.Qt.DashLine), name="Na Continuum")

            # Subtract the continuum
            self.Plot.setLabel('left', 'Intensity')
            self.Plot.setLabel('bottom', 'Wavelength (nm)')
            self.Plot.setXRange(np.min(wavelengths), np.max(wavelengths))
            self.Plot.setYRange(0, np.max(intensities) * 1.2)


            # Trying to show the lines between the top and bottom points
            # Calculate average wavelength for top and bottom points
            if top_points and bottom_points:
                # Get the wavelengths corresponding to the selected points
                # (You already have selected_top_wavelengths and selected_bottom_wavelengths in calculateAverageFeIntensity)
                # If not, you can reconstruct them here:
                all_fe_points_desc = sorted(all_fe_points, key=lambda x: x[1], reverse=True)
                all_fe_points_asc = sorted(all_fe_points, key=lambda x: x[1])
                selected_top_wavelengths = []
                selected_bottom_wavelengths = []

                # Reconstruct top wavelengths
                for wave, inten in all_fe_points_desc:
                    if any(abs(wave - sel_wave) < 1.0 for sel_wave in selected_top_wavelengths):
                        continue
                    selected_top_wavelengths.append(wave)
                    if len(selected_top_wavelengths) == len(top_points):
                        break

                # Reconstruct bottom wavelengths
                for wave, inten in all_fe_points_asc:
                    if any(abs(wave - sel_wave) < 1.0 for sel_wave in selected_bottom_wavelengths):
                        continue
                    selected_bottom_wavelengths.append(wave)
                    if len(selected_bottom_wavelengths) == len(bottom_points):
                        break

                # Calculate the average wavelength for each
                avg_top_wavelength = np.mean(selected_top_wavelengths)
                avg_bottom_wavelength = np.mean(selected_bottom_wavelengths)

                # Plot the trendline between these two points
                # self.Plot.plot(
                #     [avg_bottom_wavelength, avg_top_wavelength],
                #     [bottom_n_average, top_n_average],
                #     pen=pg.mkPen(color=(30,144,255), width=2, style=QtCore.Qt.DashLine),
                #     name="Fe Trend (Lowest to Highest)"
                # )

            print("Continuum subtraction models plotted over spectrum.")


        except Exception as e:
            print(f"Error during continuum subtraction: {e}")


    def saveSubtractionInfo(self):
        """
        Calls subtractContinuum and saves its outputs to a file.
        """

        # Calculate average intensity of Fe lines first
        fe_ranges=[(420, 455), (485, 515), (520, 550)]
        all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
        max_fe_intensity = max([inten for wave, inten in all_fe_points]) if all_fe_points else 0.0
        # print(f"Maximum Fe intensity in ranges: {max_fe_intensity:.6f}")
        
        # fe_average_value = (top_n_average + bottom_n_average) / 2
        # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0
        fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0

        # Followed by continuum subtraction
        intensity_dict, fe_average_value_new, norm_intensities = self.subtractContinuum()

        # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0 # might not need the first call?
        # fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0

        if intensity_dict is None:
            print("Subtraction failed, nothing to save.")
            return

        save_path = self.SavePath_edit.text()
        if not save_path:
            print("Save path is not set. Please choose a save path.")
            return

        file_name = os.path.join(save_path, "saved_subtraction_info_horizontal.txt")
        try:
            with open(file_name, 'w') as f:
                f.write("Subtraction Info\n")
                f.write("================\n")
                for fe_range, average_intensity in fe_averages.items():
                    f.write(f"Fe range {fe_range}: Average intensity = {average_intensity:.6f}\n")
                f.write(f"\nTop Fe intensities (with exclusion of +/- 1.0 nm):\n")
                for i, val in enumerate(top_points):
                    f.write(f"  {i+1}: {val:.6f}\n")
                f.write(f"\nAverage of top {len(top_points)} Fe intensities: {top_n_average:.6f}\n")

                f.write(f"\nBottom Fe intensities (with exclusion of +/- 1.0 nm):\n")
                for i, val in enumerate(bottom_points):
                    f.write(f"  {i+1}: {val:.6f}\n")
                f.write(f"\nAverage of bottom {len(bottom_points)} Fe intensities: {bottom_n_average:.6f}\n")
                # f.write(f"\nAverage Fe intensity: {fe_average_value:.6f}\n")
                f.write(f"\nMaximum Fe Intensity: {max_fe_intensity:.6f}\n")
                f.write(f"Intensity at 518 nm (Mg): {intensity_dict['Mg_518nm']:.6f}\n")
                f.write(f"Intensity at 589 nm (Na): {intensity_dict['Na_589nm']:.6f}\n\n")
                # f.write(f"Sum of intensities (Mg + Na + Fe): {intensity_dict['Mg_518nm'] + intensity_dict['Na_589nm'] + fe_average_value:.6f}\n\n")
                f.write(f"Sum of intensities (Mg + Na + Fe): {intensity_dict['Mg_518nm'] + intensity_dict['Na_589nm'] + max_fe_intensity:.6f}\n\n")
                f.write(f"Normalized intensities:\n")
                for k, v in norm_intensities.items():
                    f.write(f"  {k}: {v:.6f}\n")
            print(f"Subtraction info saved to {file_name}")
        except Exception as e:
            print(f"Error saving subtraction info: {e}")

    ## For the fifth tab used for background subtraction ; will keep functions only here eventually
    
    def chooseSavePath_bg(self):
        print('test')
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.Directory)

        if dlg.exec():
            SaveDirectory = str(dlg.getExistingDirectory(self, 'Select Directory'))
            print(SaveDirectory)
            self.SavePath_edit_bg.setText(SaveDirectory)


    def subtractContinuum_bg(self):
        """
        Perform linear continuum subtraction around Mg (518nm) and Na (589nm).
        """
        try:
            # Define wavelength ranges for continuum subtraction
            mg_continuum_left = (485, 515)  # Left continuum range for Mg
            mg_continuum_right = (520, 550)  # Right continuum range for Mg
            na_continuum_left = (570, 585)  # Left continuum range for Na
            na_continuum_right = (593, 608)  # Right continuum range for Na

            # Accessing user-specified regions to subtract from
            # mg_left = self.MgLeft_spin.value()
            # mg_right = self.MgRight_spin.value()
            # na_left = self.NaLeft_spin.value()
            # na_right = self.NaRight_spin.value()

            # Calculate average intensity of Fe lines
            fe_ranges=[(420, 455), (485, 515), (520, 550)]
            all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
            
            fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0
            # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0
            # fe_average_value = (top_n_average + bottom_n_average) / 2 if top_n_average and bottom_n_average else 0.0
            print(f"Top N Fe average intensity: {top_n_average:.6f}")
            print(f"Bottom N Fe average intensity: {bottom_n_average:.6f}")
            print(f"Overall Fe average intensity: {fe_average_value:.6f}")

            # Access values from the fe_averages dictionary
            for fe_range, average_intensity in fe_averages.items():
                print(f"Fe range {fe_range}: Average intensity = {average_intensity:.6f}")

            # Calculate the average of Fe intensities
            #fe_average_value = np.mean(list(fe_averages.values())) if fe_averages else 0.0
            print(f"Average Fe intensity: {fe_average_value:.6f}")

            # Extract the spectrum data
            wavelengths = np.array(self.spectrumX)
            intensities = np.array(self.spectrumY_resp) 

            # Fit linear continuum for Mg
            mg_left_mask = (wavelengths >= mg_continuum_left[0]) & (wavelengths <= mg_continuum_left[1])
            mg_right_mask = (wavelengths >= mg_continuum_right[0]) & (wavelengths <= mg_continuum_right[1])
            mg_continuum_wavelengths = np.concatenate([wavelengths[mg_left_mask], wavelengths[mg_right_mask]])
            mg_continuum_intensities = np.concatenate([intensities[mg_left_mask], intensities[mg_right_mask]])
            # # For trendline
            # mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
            # mg_baseline = np.polyval(mg_fit, wavelengths)

            # For horizontal line
            mg_level = np.median(mg_continuum_intensities)
            mg_baseline = np.full_like(wavelengths, mg_level)

            # Fit linear continuum for Na
            na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
            na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
            na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
            na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
            # # For trendline
            # na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
            # na_baseline = np.polyval(na_fit, wavelengths)

            # For horizontal line
            na_level = np.median(na_continuum_intensities)
            na_baseline = np.full_like(wavelengths, na_level)


            # Subtract the continuum
            corrected_intensities = intensities.copy()
            corrected_intensities[(wavelengths >= 485) & (wavelengths <= 550)] -= mg_baseline[(wavelengths >= 485) & (wavelengths <= 550)]
            corrected_intensities[(wavelengths >= 570) & (wavelengths <= 608)] -= na_baseline[(wavelengths >= 570) & (wavelengths <= 608)]

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

            norm_intensities['Fe'] = round(fe_average_value / intensity_sum, 6) if intensity_sum != 0 else 0.0
            norm_intensities['Mg'] = round(intensity_dict['Mg_518nm'] / intensity_sum, 6) if intensity_sum != 0 else 0.0
            norm_intensities['Na'] = round(intensity_dict['Na_589nm'] / intensity_sum, 6) if intensity_sum != 0 else 0.0
            
            print(f"Normalized intensities: {norm_intensities}")
            # Store the normalized intensities in the spectral object
            self.spectral.normalized_intensities = norm_intensities

            # At the end, instead of just printing, return the results:
            return intensity_dict, fe_average_value, norm_intensities

        except Exception as e:
            print(f"Error during continuum subtraction: {e}")
            return None, 0.0, {}
        
    def saveSubtractionInfo_bg(self):
        """
        Calls subtractContinuum and saves its outputs to a file.
        """

        # Calculate average intensity of Fe lines first
        fe_ranges=[(420, 455), (485, 515), (520, 550)]
        all_fe_points, fe_averages, top_points, top_n_average, bottom_points, bottom_n_average = self.calculateAverageFeIntensity(fe_ranges)
        # fe_average_value = (top_n_average + bottom_n_average) / 2
        # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0
        fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0

        # Followed by continuum subtraction
        intensity_dict, fe_average_value_new, norm_intensities = self.subtractContinuum()

        # fe_average_value = np.mean(all_fe_points) if len(all_fe_points) > 0 else 0.0 # might not need the first call?
        # fe_average_value = np.mean([inten for wave, inten in all_fe_points]) if len(all_fe_points) > 0 else 0.0

        if intensity_dict is None:
            print("Subtraction failed, nothing to save.")
            return

        save_path = self.SavePath_edit_bg.text()
        if not save_path:
            print("Save path is not set. Please choose a save path.")
            return
    
        

        file_name = os.path.join(save_path, "saved_subtraction_info_horizontal.txt")
        try:
            with open(file_name, 'w') as f:
                f.write("Subtraction Info\n")
                f.write("================\n")
                for fe_range, average_intensity in fe_averages.items():
                    f.write(f"Fe range {fe_range}: Average intensity = {average_intensity:.6f}\n")
                f.write(f"\nTop Fe intensities (with exclusion of +/- 1.0 nm):\n")
                for i, val in enumerate(top_points):
                    f.write(f"  {i+1}: {val:.6f}\n")
                f.write(f"\nAverage of top {len(top_points)} Fe intensities: {top_n_average:.6f}\n")

                f.write(f"\nBottom Fe intensities (with exclusion of +/- 1.0 nm):\n")
                for i, val in enumerate(bottom_points):
                    f.write(f"  {i+1}: {val:.6f}\n")
                f.write(f"\nAverage of bottom {len(bottom_points)} Fe intensities: {bottom_n_average:.6f}\n")
                f.write(f"\nAverage Fe intensity: {fe_average_value:.6f}\n")
                f.write(f"Intensity at 518 nm (Mg): {intensity_dict['Mg_518nm']:.6f}\n")
                f.write(f"Intensity at 589 nm (Na): {intensity_dict['Na_589nm']:.6f}\n\n")
                f.write(f"Sum of intensities (Mg + Na + Fe): {intensity_dict['Mg_518nm'] + intensity_dict['Na_589nm'] + fe_average_value:.6f}\n\n")
                f.write(f"Normalized intensities:\n")
                for k, v in norm_intensities.items():
                    f.write(f"  {k}: {v:.6f}\n")
            print(f"Subtraction info saved to {file_name}")
        except Exception as e:
            print(f"Error saving subtraction info: {e}")

# Added to the savePlot function

# --- Add continuum model lines if the toggle is checked ---
        if hasattr(self, "ShowContinuumModel_check") and self.ShowContinuumModel_check.isChecked():

            try:
                mg_continuum_left = (510, 515)
                mg_continuum_right = (520, 525)
                na_continuum_left = (580, 585)
                na_continuum_right = (593, 598)

                wavelengths = np.array(self.spectrumX)
                intensities = np.array(self.spectrumY_resp)

                # Mg continuum
                mg_left_mask = (wavelengths >= mg_continuum_left[0]) & (wavelengths <= mg_continuum_left[1])
                mg_right_mask = (wavelengths >= mg_continuum_right[0]) & (wavelengths <= mg_continuum_right[1])
                mg_continuum_wavelengths = np.concatenate([wavelengths[mg_left_mask], wavelengths[mg_right_mask]])
                mg_continuum_intensities = np.concatenate([intensities[mg_left_mask], intensities[mg_right_mask]])
                mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
                mg_baseline = np.polyval(mg_fit, wavelengths)
                mg_region_mask = (wavelengths >= 510) & (wavelengths <= 525)
                ax.plot(wavelengths[mg_region_mask], mg_baseline[mg_region_mask], 
                        color='red', linestyle='--', linewidth=2, label='Mg Background Continuum')

                # Na continuum
                na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
                na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
                na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
                na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
                na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
                na_baseline = np.polyval(na_fit, wavelengths)
                na_region_mask = (wavelengths >= 580) & (wavelengths <= 598)
                ax.plot(wavelengths[na_region_mask], na_baseline[na_region_mask], 
                        color='black', linestyle='--', linewidth=2, label='Na Background Continuum')
            except Exception as e:
                print(f"Could not plot continuum models in savePlot: {e}")

        print("ShowContinuumModel_check:", self.ShowContinuumModel_check.isChecked())


# Added to the plotMeasuredSpec function

# # After plotting the main spectrum, add the continuum models:
        # This shows it in the plot, but it is not used in the final plot
        # # Only plot continuum models if the toggle is checked

        if hasattr(self, "ShowContinuumModel_check") and self.ShowContinuumModel_check.isChecked():
        
            try:
                # Define continuum regions
                mg_continuum_left = (510, 515)
                mg_continuum_right = (520, 525)
                na_continuum_left = (580, 585)
                na_continuum_right = (593, 598)

                wavelengths = np.array(self.spectrumX)
                intensities = np.array(self.spectrumY_resp)

                # Mg continuum
                mg_left_mask = (wavelengths >= mg_continuum_left[0]) & (wavelengths <= mg_continuum_left[1])
                mg_right_mask = (wavelengths >= mg_continuum_right[0]) & (wavelengths <= mg_continuum_right[1])
                mg_continuum_wavelengths = np.concatenate([wavelengths[mg_left_mask], wavelengths[mg_right_mask]])
                mg_continuum_intensities = np.concatenate([intensities[mg_left_mask], intensities[mg_right_mask]])
                mg_fit = np.polyfit(mg_continuum_wavelengths, mg_continuum_intensities, deg=1)
                mg_baseline = np.polyval(mg_fit, wavelengths)
                mg_region_mask = (wavelengths >= 510) & (wavelengths <= 525)
                self.Plot.plot(wavelengths[mg_region_mask], mg_baseline[mg_region_mask],
                            pen=pg.mkPen(color=(255, 0, 0), width=2, style=QtCore.Qt.DashLine), name="Mg Background Continuum")

                # Na continuum
                na_left_mask = (wavelengths >= na_continuum_left[0]) & (wavelengths <= na_continuum_left[1])
                na_right_mask = (wavelengths >= na_continuum_right[0]) & (wavelengths <= na_continuum_right[1])
                na_continuum_wavelengths = np.concatenate([wavelengths[na_left_mask], wavelengths[na_right_mask]])
                na_continuum_intensities = np.concatenate([intensities[na_left_mask], intensities[na_right_mask]])
                na_fit = np.polyfit(na_continuum_wavelengths, na_continuum_intensities, deg=1)
                na_baseline = np.polyval(na_fit, wavelengths)
                na_region_mask = (wavelengths >= 580) & (wavelengths <= 598)
                self.Plot.plot(wavelengths[na_region_mask], na_baseline[na_region_mask],
                            pen=pg.mkPen(color=(0, 0, 0), width=2, style=QtCore.Qt.DashLine), name="Na Background SContinuum")
                
            except Exception as e:
                print(f"Could not plot continuum models: {e}")
        



# Put these into the Ui Class under the init function

 # Save fitted elements
        self.SaveFittedElements_button.clicked.connect(self.saveFittedElements)
        #MJM
        self.ChooseSavePath_button.clicked.connect(self.chooseSavePath)

        # Fe Average Button
        self.AverageFe_button.clicked.connect(self.calculateAverageFeIntensity)

        # Continuum subtraction
        self.ContinuumSubtract_button.clicked.connect(self.subtractContinuum)

        # Toggle for showing continuum models
        self.ShowContinuumModel_check.toggled.connect(self.refreshPlot)

        # Save continuum subtraction info
        self.SaveSubtractionInfo_button.clicked.connect(self.saveSubtractionInfo)

        # For subtraction tab
        self.ChooseSavePath_button_bg.clicked.connect(self.chooseSavePath_bg)

        # For subtraction tab
        self.SaveSubtractionInfo_button_bg.clicked.connect(self.saveSubtractionInfo_bg)

        # For subtraction tab
        self.ContinuumSubtract_button_bg.clicked.connect(self.subtractContinuum_bg)

        # Continuum model overlay
        self.ContinuumModel_button.clicked.connect(self.showContinuumModel)

# initialize the GuralSpectral object
        self.spectral = spectral_library.GuralSpectral(10000, 4500, None, None, None, None)      


    # Need to add buttons for the above Ui functions into the UI as well


# At the top of file where modules are imported

# ...existing code...
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# ...existing code...

from PyQt5.QtWidgets import QCheckBox  # Already imported, but for clarity
