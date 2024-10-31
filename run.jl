# run.jl - This script is used to run the scripts in the project. Set the raw_data_dummy to
# true if you want to run the scripts to produce the analysis data from the raw_data. If 
# choosing this option you must take steps to download and format the raw data as outlined
# in the ReadMe.pdf file in the root of the project.

############################################################################################
#### Recommended Requirements ##############################################################
############################################################################################
# Memory: ≥ 32 GBs of RAM. Hard Drive Space: ≥ 20 GBs of free space. Processor: ≥ 4 cores.
# The code is written to be OS agnostic with file paths. Install Julia 1.10.4 and set the 
# path. In Visual Studio Code, hit ctrl + , to open the settings panel. Search for 
# "Julia: executable path" and set to the path of the Julia executable. On Windows, it will 
# look something like "C:\Users\davle\AppData\Local\Programs\Julia-1.10.4\bin\julia.exe".
# Also search for "Julia: Num Threads". Click "edit in settings.json".
# Edit the following line to the number of cores you want to use, ""julia.NumThreads": 6,",
# you would change "6" to the number of cores you want to use.  If you plan producing the
# map figure, you must install R version 4.3.1 (2023-06-16 ucrt) -- "Beagle Scouts" and
# include it in the system path. The code will install all packages neccesary. Lastly, if 
# setting use_rcall to true, which will run the R script to produce the map, you must set
# the "R_HOME" environment variable to the path of the R installation. See
# https://juliainterop.github.io/RCall.jl/dev/installation/#Customizing-the-R-installation-using-R_HOME
# for further details.
############################################################################################

############################################################################################
#### Bechmarking Specifications ############################################################
############################################################################################
# The code was run on a desktop pc equipped with an AMD Ryzen 5 56700 6-core processor,
# 64 GB of ram, and running Windows 11 Pro. This script was sequentially run using Julia 
# 1.10.4 using Visual Studio Code. The code was run with the following parameters:
# produce_raw_data = true, save_intermediate_files = false, produce_map = true, and 
# use_rcall = true.
############################################################################################
##### Benchmarking Results #################################################################
############################################################################################
# 1 - saving_intermediate_files = false 
#   - The code completed in 2655.04 seconds ≈ 44.25
#     The size of the directory created was 8.31gb and compressed to 906mb.
#
# 2 - saving_intermediate_files = true
#   - The code completed in 3333.16 seconds ≈ 55.55 minutes
#     The size of the directory created was 28.2gb and compressed to 1.81gb
#
############################################################################################
#### Code Parameters ##############################################################
############################################################################################
#If you are not in the root directory, change path to the root directory
#cd(joinpath(path, to, root))
seed_num = 12431413412

produce_raw_data = true # Set to true to run the scripts to produce the analysis data from the
                      # raw_data
                    
save_intermediate_files = false #Set to true to write intermediate files. Adds 13gb.

produce_map = true # Set to true to produce the map of the samples. If use_rcall is set to
                   # true, the R script will be run to produce the map. If use_rcall is set
                   # to false, then the file processed_data/df_map.csv will be written and 
                   # the user can use the R script, create_sample_map.R to produce the map.

use_rcall = true # Set to true to use RCall to run the R script to produce the map.

############################################################################################
#### Set-up the project encironment ########################################################
############################################################################################
using Pkg
# Activate the current environment
Pkg.activate("./code/Bunching")
# Install and instantiate the project
Pkg.instantiate()
# Precompile the packages
Pkg.precompile() # Precompilation time: 242s ≈ 4m; Size on Disk: 630mb
using Bunching

@time begin "Code completed in: "
   
############################################################################################
#Bench mark entire code w/ produce_raw_data = true, save_intermediate_files = true
# produce_map = true, use_rcall = true. 
############################################################################################
#### Produce analysis data #################################################################
############################################################################################
if produce_raw_data
    #Benchmarking individually with save_intermediate_files = true
    df = clean_raw_pluto_data(save_intermediate_files) #time: 1738.93s ≈ 30m
    df = transform_cleaned_pluto_data(df, save_intermediate_files) #time: 1473.27s ≈ 25m
    generate_analysis_subsamples(df)  # time: 3.3s
end
############################################################################################

############################################################################################
#### Generate tables and figures ##########################################################
############################################################################################
generate_tables_figures(seed_num)
if produce_map
    create_sample_map(use_rcall)
end
############################################################################################
end