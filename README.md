# ciceroscm
Python version of the CICERO-SCM simple climate model/emulator

## Running
To run the model, copy <code>run_full_forcing.py</code> or <code>run_scm.py</code>
for a forcing, or full emissions run from the <code>scripts</code> directory. The forcing run currently
supports only files with single total forcing series supplemented
with solar and volcanic forcing time series.
There is currently no support for pure concentration runs, only full
emissions.
The scripts show you how to specify input data and output data placement.
The output data will be placed in a folder called output_test in the
folder from which you run the script.
More input data on the appropriate format (at least for emission runs) can be found in <code>/div/amoc/CSCM/SCM_Linux_v2019/RCMIP/input/</code> on amoc for internal use.

## Instance configurations
When a new instance of the CICERO-SCM class is created the dictionary cfg needs to be sent as a parameter, detailing the configurations of the instance. All options that end in _file can be exchanged with the same parameter ending in _data to send in data directly. See details on the data format in the input_handler module. Configuration options are:
* nystart - the start year of the run
* nyend - the end year of the run
* emstart - the year to start the run with emissions
* idtm - optional parameter to tune the number of subyearly steps in the concentrations_emissions_handler. Default is 24. Should probably not be the first parameter you want to start playing with.
* sunvolc - an optional parameter to include solar an volcanic forcing. If included and equal to 1 a set of such forcing series will be included.
* rf_sun_file - optional path to file with solar data, only read if sunvolc is 1. If sunvolc is 1 and this parameter is not set, a default file will be used.
* rf_volc_file - optional path to file with hemispherically symmetric volcanic forcing data. If you prefer, you can send rf_volc_n_file and rf_volc_s_file for separate data for each hemisphere, but then a global file must not be sent, as it will override the hemispherically split files when present. The volcanic data can be on columns, with monthly data, on one yearly column, or on any other periodic split per column per year (i.e. seasonal, half yearly, every four months). If sunvolc is not 1, all of these will be ignored. If sunvolc is 1 and none and no volcanic forcing data is indicated by the user, a default file will be used.
* gaspam_file - path to file of gases to include with units, forcing factors etc (mandatory).
* concentrations_file - path to file with concentrations time series (mandatory if not forcing run)
* emissions_file - path to file with emissions time series (mandatory if not forcing run)
* nat_ch4_file- optional path to file where natural emissions timeseries for methane can be found. Default file will be used if not given
* nat_n2o_file- optional path to file where natural emissions timeseries for n2o can be found. Default file will be used if not given
* forc_file - path to file with forcing time series, if this is sent the run will be a forcing run, and none of the emission and concentration related options will be relevant. The file can be a single column of numbers of total forcing, it will be assumed to run from whatever startyear you set, or a comma separated file, with 'year' as first column, followed by either hemispherically split forcing under headings "FORC_NH" and "FORC_SH", or columns per various forcing components. (At the moment you cannot include hemispherical split along with several components)
* conc_run - Set this to True and have a concentration driven run. You will still need to provide an emission file, as some species forcings (such as ozone) are calculated from emissions after emstart.
* perturb_em_file - path to file with emission perturbations to be added to the emissions from the emissions file, the format for this file is shown in the file in test/test_data/pertem_test.txt
* perturb_forc_file - path to file with forcings to be added after forcings from emissions and concentrations have been calculated, the format for this file is shown in the file in test/test_data/pertforc_test.txt

## Options for run
With a CICEROSCM instance in place, you are ready to start runs with various parameter configurations, using the input files as set by the instance configuration

### Run configurations
* output_folder - name of or path of file wher output from the run is stored (at the moment this will always be assumed to be laying under the directory from which the code is run)
* output_prefix - prefix to output filenames
* make_plot - if set to True plots of the output are made and saved to a subfolder in the output_folder.

### Parameter sets
Physical parameters to the model is divided in two parametersets each of which are sent as two seperate dictionaries to the run call.
* One parameterset pamset_udm for the upwelling diffusion model
* One parameterset pamset_emiconc for emissions and concentrations
If the parametersets are not provided, a default parameterset is used
If one or more parameters are not provided as part of the parameterset, these parameters will be set to the default values
#### pamset_udm
The upwelling diffusion model (which is needed for all runs) takes the following parameters.(Default value in paranthesis):
* rlamdo (16.0) - Air-sea heat exchange parameter <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2\mathrm{K}}">, valid range 5-25
* akapa (0.634) - Vertical heat diffusicity <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{cm}^2}{\mathrm{s}}">, valid range 0.06-0.8
* cpi (0.4) - Polar parameter, temperature change ratio polar to nonpolar region, unitless, valid range 0.161-0.569
* W (4.0) - Vertical velocity, upwelling rate <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{m}}{\mathrm{yr}}">, valid range 0.55-2.55
* beto (3.5) - Ocean interhemispheric heat exchange coefficient <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2\mathrm{K}}">, valid range 0-7
* threstemp (7.0) - Scales vertical velocity as a function of mixed layer temperature, unitless. Set to 0 if you don't want to include this parameter.
* lambda (0.540) - Equilibrium climate sensitivity divided by 2xCO2 radiative forcing <img src="https://render.githubusercontent.com/render/math?math=\large \left( 3.71 \frac{\mathrm{W}}{\mathrm{m}^2} \right)">
* mixed (60.0) - Mixed layer depth, m, valid range 25-125
* foan (0.61) - Fraction of Northern hemisphere covered by ocean
* foas (0.81) - Fraction of Southern hemisphere covered by ocean
* ebbeta (0.0) - Atmospheric interhemispheric heat exchange (not currently used)
* fnso (0.7531) - Ratio between ocean areas in Northern and Southern hemispheres
* lm (40) - Number of vertical layers (below the mixed layer each layer is at each 100 m depth)
* ldtime (12) - Number of subyearly timesteps

#### pamset_emiconc
The concentration and emission parameterset (which is needed for emission runs) takes the following parameters. (Default value in paranthesis):

* qbmb (0.03) - Biomass burning aerosol RF in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qo3 (0.4) - Tropospheric ozone RF in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qdirso2 (-0.457) - Direct RF sulphate in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qindso2 (-0.514) - Indirect RF sulphate in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qbc (0.200) - BC (fossil fuel + biofuel) RF in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qoc (-0.103) - OC (fossil fuel + biofuel) RF in ref_yr, <img src="https://render.githubusercontent.com/render/math?math=\large \frac{\mathrm{W}}{\mathrm{m}^2}">
* qh2o_ch4 (0.091915) - Stratospheric water vapour ERF ratio to methane ERF
* ref_yr (2010) - Reference year for the above forcing values. To construct radiative forcing time series, these forcing values are scaled using emssions. The forcing in the reference year is equal to the forcing value set by the above parameters
* idtm (24) - Number of subyearly timesteps for calculation of CO2 concentrations from emissions.

## Example scripts
The scripts folder contains various example scripts that can be used to see how to set up various types of runs. The start of all of them adds the necessary parts for the file to run with the module. If you want to run from somewhere else you will need to edit the <code>sys.path.append</code> command so it points to where you've stored the src directory of this repository.
* <code>run_scm.py</code> runs a simple emissions run with ssp245 data from 1900 to 2050
* <code>run_full_forcing.py</code> runs a 1 percent CO2 increase forcing with default solar and volcanic forcing from 1750 to 2100
* <code>run_perturbations.py</code> shows runs like that in <code>run_scm.py</code> with emissions and forcing perturbations
* <code>run_full_emissions_profile.py</code> runs an ssp245 emissions run from 1750 to 2100 with a profiler, so you can see what parts of the code is more time consuming
* <code>run_full_forcing_profile.py</code> is like <code>run_full_emissions_profile.py</code> but for a pure forcing run
* <code>run_full_change_all_pams.py</code> is an emissions ssp245 run from 1750 to 2100 which shows how to set all the parameters for both the upwelling diffusion model and for the concentrations emissions handler.
* <code>run_ssps_local.py</code> runs through all scenarios on on amoc, this script will only work on amoc or qbo, but can show how to loop through elsewhere, just remember to change paths.

## Jupyter Notebooks
The notebooks folder provides simple working examples to run the model within a Jupyter environment, and plot example output.  Installation instructions for installing Jupyterlab can be found at https://jupyter.org/install
<code>CSCM_example_textinput.ipynb</code> runs a simple emissions run with ssp245 data from 1900 to 2050, using input data text files stored in the <code>tests/test-data</code> folder
<code>CSCM_example_directinput.ipynb</code> illustrates an interactive case, where ssp data is read into the environment and passed directly to the model

### prescripts
Inside the scripts folder is a folder called prescripts. It contains scripts that show how to prepare perturbation files for a run and two example datafiles. And includes scripts to prepare natural emissions files.

## Development
* To start developing make sure you have a github account and that you are part of the ciceroOslo team.
* If you haven't already, [setup your github account with an ssh key](https://docs.github.com/en/enterprise-server@3.0/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
* Find a suitable place where you want to start your developement either on windows or under /div/nobakcup/users/yourusername/ and change to that directory
* Once in the preferred directory with a terminal do:
> <code>git clone git@github.com:ciceroOslo/ciceroscm.git</code>
* To make your own branch (which you really should)
> <code>git checkout -b your-cool-branch-name</code>
* Whenever you log in or want to check stuff
> <code>git status</code>

It will tell you the branch you are on, changes since last etc
* To commit your code changes
> <code>git add path-of-file-that-changed</code>

Repeat this for all the files that you would want to commit the changes for
> <code>git commit -m "A small message to describe the changes"</code>

> <code>git push</code>

(The last one is to push the changes to the github version. The first time youi do this on a new branch you will need to set where to push to, but how to do that will be suggested when you just do git push)
* To get new changes that have happened on the main branch is always good before you commit. To do so do:
> <code>git checkout main</code>

> <code>git pull</code>

> <code>git checkout your-cool-branch-name</code>

> <code>git merge main</code>

If all goes well this will fill your terminal with a merge message in your default editor, which is likely vim. The message there is likely ok as it is, so to just use that as a commit message for the merge type: <code>:wq</code> which will just save and quit vim and complete the merge with the original commit message.

Then finally just push your code to the web.

> <code>git push</code>

The last part is just to pushed this new version of your branch again

### Test suite and environment
The code comes with a suite of tests and tools. To use this you must do:
> <code>make first-venv</code>

> <code>make virtual-environment</code>

This should only be necessary the first time you setup the code
You can load this environment with
> <code>source venv/bin/activate</code>

Later to update you should do:
> <code>make virtual-environment</code>

Or if you know you need updates, but aren't getting them:
> <code>make clean</code>

> <code>make virtual-environment</code>

After this you should be able to run the automatic tests
> <code>make test</code>

Will only run the tests
> <code>make checks</code>

Will run the tests and formatting tests

Before your code branch can be merged into the main code, it has to pass all the tests
(The makefile also has an option to run only the formatting checks)
Tests are located in tests in tests/test-data/ data for testing against fortran runs and test input data are stored. In tests/unit there are unit tests for certain methods. In test/integration there are integration tests of the code, comparing the results to fortran.
When you develop new code, try to think about what can be done to test and validate that your code does what you expect it to do, and try to integrate such tests into the automatic testing scheme.

## General code flow
The code consists of four modules
* ciceroscm takes in an sorts inputs, is what gets called, and loops over the years and calls the other methods. It also outputs temperature and ocean data.
* upwelling_diffusion_method is the energy budgeting method that takes forcing to temperature, ocean heat content etc. It gets called and delivers results to ciceroscm.
* concentration_emissions_handler takes care of calculating its way from emissions to concentrations to forcing. It gets called every year, but saves it's results internally and only returns the forcing. It also has an output method of it's own to produce the emission, concentration and forcing files from the run
* _utils is just a method to put common utilities in. At the moment it has only one method that can check whehter a parameterset includes the expected values and putting in default values if not.
* perturbations.py handles and adds perturbations to either forcing or emissions per species.
* make_plots makes plots if plotting functionality is invoked.
* input_handler takes care of reading in files or data, and has various file reading methods.
