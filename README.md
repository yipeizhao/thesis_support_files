# thesis_support_files
Including all files that are used to build my master thesis
Files/folders descriptions:
* bus: Includes all bus networks(raw/processed/modified)
* env: A python virtual environment that allows you to run the exactly version of python and python packages as mine. Detailed description: https://docs.python.org/3/library/venv.html
* figures: Generated figures
* real_networks: Includes all real networks(raw/processed)
* result: Includes the results of rewiring on real networks
* Beijing_bus_process.py: Transfer the raw Beijing bus network to a usable network
* Complexity.py: Includes all complexity meaures
* configuration_model.r: Used to monitor the behaviour of gamma
* figures_generator.ipynb: Used to generate most figures
* modify_bus.py: Used to modify bus networks
* real_networks_preprocess.ipynb: Used to convert raw bus networks(other networks need to be in the form stated in the data section in the thesis, exclude GBPT_train, and Beijing, London bus network)
* renv.lock: Allows you to run exactly same version of R and R packages as mine. Detailed description: https://rstudio.github.io/renv/articles/renv.html
* rewiring.ipynb: Used to apply the rewiring techniques
* R_functions.py: Used to supply configuration_model.r
* utilities.py: Includes various functions to support the project
* version.ipynb: Includes the version of python and the used python packages 
