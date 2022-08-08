ciceroscm
==========

A python version of the cicero-scm simple climate model

License
-------

.. sec-begin-license
Apache-2.0 license

.. sec-end-license
.. sec-begin-long-description

Background
----------
This is the python version of the CICERO-SCM simple climate model/emulator.

The Cicero Simple Climate model was first developed over 20 years ago and has been in use at Cicero with various updated ever since. Broadly speaking, it goes from global emissions to concentrations, to radiative forcing to global temperature change. The model can beshort circuited to start at any of these inputs, running directly from forcing or concentrations rather than all the way. However, for shortlived forcers, forcing is calculated directly from emisssions, so emissions need to be supplied even in a concentration driven run.

The model includes a carbon cycle model introduced in a 1996 paper by Joos et al. The interplay between CO2, CH4 and N2O has been updated according to the 2016 work of Etminan et al. Otherwise, simplified impulse response functions take Tg emissions to concentrasions using lifetimes.

Natural emissions of CH4 and N2O should be time variable and tuned in order to fit historical evolution.

The forcing to temperature calculation is acheived via the 1992 Schlesinger upwelling diffusion model, calculating energy balance in a default set of 40 layered ocean with a box atmosphere for each hemisphere. Both forcing calculations and energy budget can be tuned using various user supplied parameters.

Currently, this is a fairly faithful python rendering of the pre-existing fortran model, albeit with some added flexibility, some AR6 updates, and a couple of bugfixes.



.. sec-end-long-description

.. sec-begin-installation

Installation
------------
CiceroSCM can not be installed via pip or conda yet.

.. sec-end-installation

Documentation
-------------

Documentation in docs folder not properly setup elsewhere yet.

Contributing
------------

Please see the developement section of docs. At the moment only internal contributers.

.. sec-begin-links

.. sec-end-links
