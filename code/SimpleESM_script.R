### Calculation of SOC (and N) stocks, in t/ha, at fixed depth (FD) and Equivalent Soil Mass (ESM) ###

# Authors: Fabien FERCHAUD and Florent CHLEBOWSKI
# Contact: Fabien FERCHAUD (INRAE) - fabien.ferchaud@inrae.fr
# Date: 2022-12-16

# Packages -------
library(here)

# Options -----------------------------------------------------------------

# Name of the input .xlsx file
input_file_name <- here("data_raw", "Test_SimpleESM_kansas_max.xlsx")

# Name of the output directory
output_directory_name <- here("data_processed")

# Option for the reference soil mass ("manual" or "auto")
RefM_option <- "manual"

# Option for calculations - elements: one or two elements ("SOC_only" or "SOC_and_N")
E_calc_option <- "SOC_only"

# Option for calculations - isotopes: 13C, or 13C and 15N, or not ("13C" or "13C_15N" or "no")
I_calc_option <- "no"


# Call SimpleESM function -------------------------------------------------

source(here("code","SimpleESM_function.R"))

SimpleESM(input_file_name, output_directory_name, RefM_option, E_calc_option, I_calc_option)
