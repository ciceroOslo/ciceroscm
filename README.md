# ciceroscm
Python version of the CICERO-SCM simple climate model/emulator

## Running
To run the model, copy the script run_full_forcing.py or run_scm.py
for a forcing, or full emissions run. The forcing run currently
supports only files with single total forcing series supplemented
with solar and volcanic forcing time series.
There is currently no support for pure concentration runs, only full
emissions.
The scripts show you how to specify input data and output data placement.
The output data will be placed in a folder called output_test in the
folder from which you run the script.
More input data on the appropriate format (at least for emission runs) can be found in /div/amoc/CSCM/SCM_Linux_v2019/RCMIP/input/ on amoc for internal use.

## Paramters
Physical parameters to the model is divided in two parametersets
* One parameterset pamset_udm for the upwelling diffusion model
* One parameterset pamset_emiconc for emissions and concentrations
If the parametersets are not provided, a default parameterset is used
If one or more parameters are not provided as part of the parameterset, these parameters will be set to the default values
The upwelling diffusion model (which is needed for all runs) takes the following parameters.(Default value in paranthesis):
* rlamdo (16.0)
* akapa (0.634)
* cpi (0.4)
* W (4.0)
* beto (3.5)
* threstemp (7.0)
* lambda (0.540)
* mixed (60.0)
* foan (0.61)
* foas (0.81)
* ebbeta (0.0)
* fnso (0.7531)
* lm (40)
* ldtime (12)
The concentration and emission parameterset (which is needed for emission runs) takes the following parameters. (Default value in paranthesis):
* lamb (0.8)
* qbmb (0.03)
* qo3 (0.4)
* qdirso2 (-0.457)
* qindso2 (-0.514)
* qbc (0.200)
* qoc (-0.103)
* ref_yr (2010)
* idtm (24)

## Developement
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
* upwelling_diffusion_method is the energy budgeting method that takes forcing to temperature, ocean heat content sea level rise etc. It gets called and delivers results to ciceroscm.
* concentration_emissions_handler takes care of calculating its way from emissions to concentrations to forcing. It gets called every year, but saves it's results internally and only returns the forcing. It also has an output method of it's own to produce the emission, concentration and forcing files from the run
* _utils is just a method to put common utilities in. At the moment it has only one method that can check whehter a parameterset includes the expected values and putting in default values if not.
