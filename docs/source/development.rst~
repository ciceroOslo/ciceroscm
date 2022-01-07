.. development:

Development
===========

If you're interested in contributing to OpenSCM-Runner, we'd love to have you on board!
This section of the docs will (once we've written it) detail how to get setup to contribute and how best to communicate.

.. contents:: :local:

Contributing
------------

All contributions are welcome, some possible suggestions include:

- tutorials (or support questions which, once solved, result in a new tutorial :D)
- blog posts
- improving the documentation
- bug reports
- feature requests
- pull requests

Please report issues or discuss feature requests in the `OpenSCM-Runner issue tracker`_.
If your issue is a feature request or a bug, please use the templates available, otherwise, simply open a normal issue :)

As a contributor, please follow a couple of conventions:

- Create issues in the `OpenSCM-Runner issue tracker`_ for changes and enhancements, this ensures that everyone in the community has a chance to comment
- Be welcoming to newcomers and encourage diverse new contributors from all backgrounds: see the `Python Community Code of Conduct <https://www.python.org/psf/codeofconduct/>`_
- Only push to your own branches, this allows people to force push to their own branches as they need without fear or causing others headaches
- Start all pull requests as draft pull requests and only mark them as ready for review once they've been rebased onto master, this makes it much simpler for reviewers
- Try and make lots of small pull requests, this makes it easier for reviewers and faster for everyone as review time grows exponentially with the number of lines in a pull request


Getting setup
-------------

To get setup as a developer, we recommend the following steps (if any of these tools are unfamiliar, please see the resources we recommend in `Development tools`_):

#. Install conda and make
#. Run ``make virtual-environment``, if that fails you can try doing it manually

    #. Change your current directory to OpenSCM-Runner's root directory (i.e. the one which contains ``README.rst``), ``cd openscm-runner``
    #. Create a virtual environment to use with OpenSCM-Runner ``python3 -m venv venv``
    #. Activate your virtual environment ``source ./venv/bin/activate``
    #. Upgrade pip ``pip intall --upgrade pip``
    #. Install the development dependencies (very important, make sure your virtual environment is active before doing this) ``pip install -e .[dev]``

#. Make sure the tests pass by running ``make check``, if that fails the commands can be read out of the ``Makefile``


Getting help
~~~~~~~~~~~~

Whilst developing, unexpected things can go wrong (that's why it's called 'developing', if we knew what we were doing, it would already be 'developed').
Normally, the fastest way to solve an issue is to contact us via the `issue tracker <https://github.com/openscm/openscm-runner/issues>`_.
The other option is to debug yourself.
For this purpose, we provide a list of the tools we use during our development as starting points for your search to find what has gone wrong.

Development tools
+++++++++++++++++

This list of development tools is what we rely on to develop OpenSCM-Runner reliably and reproducibly.
It gives you a few starting points in case things do go inexplicably wrong and you want to work out why.
We include links with each of these tools to starting points that we think are useful, in case you want to learn more.

- `Git <http://swcarpentry.github.io/git-novice/>`_

- `Make <https://swcarpentry.github.io/make-novice/>`_

- `Conda virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_

- `Pip and pip virtual environments <https://www.dabapps.com/blog/introduction-to-pip-and-virtualenv-python/>`_

- `Tests <https://semaphoreci.com/community/tutorials/testing-python-applications-with-pytest>`_

    - we use a blend of `pytest <https://docs.pytest.org/en/latest/>`_ and the inbuilt Python testing capabilities for our tests so checkout what we've already done in ``tests`` to get a feel for how it works

- `Continuous integration (CI) <https://help.github.com/en/actions>`_ (also `brief intro blog post <https://gabrieltanner.org/blog/an-introduction-to-github-actions>`_ and a `longer read <https://dev.to/bnb/an-unintentionally-comprehensive-introduction-to-github-actions-ci-blm>`_)

    - we use GitHub CI for our CI but there are a number of good providers

- `Jupyter Notebooks <https://medium.com/codingthesmartway-com-blog/getting-started-with-jupyter-notebook-for-python-4e7082bd5d46>`_

    - Jupyter is automatically included in your virtual environment if you follow our `Getting setup`_ instructions

- Sphinx_


Other tools
+++++++++++

We also use some other tools which aren't necessarily the most familiar.
Here we provide a list of these along with useful resources.

.. _regular-expressions:

- `Regular expressions <https://www.oreilly.com/ideas/an-introduction-to-regular-expressions>`_

    - we use `regex101.com <regex101.com>`_ to help us write and check our regular expressions, make sure the language is set to Python to make your life easy!

Formatting
----------

To help us focus on what the code does, not how it looks, we use a couple of automatic formatting tools.
These automatically format the code for us and tell use where the errors are.
To use them, after setting yourself up (see `Getting setup`_), simply run ``make format``.
Note that ``make format`` can only be run if you have committed all your work i.e. your working directory is 'clean'.
This restriction is made to ensure that you don't format code without being able to undo it, just in case something goes wrong.


Buiding the docs
----------------

After setting yourself up (see `Getting setup`_), building the docs is as simple as running ``make docs`` (note, run ``make -B docs`` to force the docs to rebuild and ignore make when it says '... index.html is up to date').
This will build the docs for you.
You can preview them by opening ``docs/build/html/index.html`` in a browser.

For documentation we use Sphinx_.
To get ourselves started with Sphinx, we started with `this example <https://pythonhosted.org/an_example_pypi_project/sphinx.html>`_ then used `Sphinx's getting started guide <http://www.sphinx-doc.org/en/master/usage/quickstart.html>`_.


Gotchas
~~~~~~~

To get Sphinx to generate pdfs (rarely worth the hassle), you require `Latexmk <https://mg.readthedocs.io/latexmk.html>`_.
On a Mac this can be installed with ``sudo tlmgr install latexmk``.
You will most likely also need to install some other packages (if you don't have the full distribution).
You can check which package contains any missing files with ``tlmgr search --global --file [filename]``.
You can then install the packages with ``sudo tlmgr install [package]``.


Docstring style
~~~~~~~~~~~~~~~

For our docstrings we use numpy style docstrings.
For more information on these, `here is the full guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_ and `the quick reference we also use <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.


Releasing
---------

First step
~~~~~~~~~~

#. Test installation with dependencies ``make test-install``
#. Update ``CHANGELOG.rst``

    - add a header for the new version between ``master`` and the latest bullet point
    - this should leave the section underneath the master header empty

#. ``git add .``
#. ``git commit -m "Prepare for release of vX.Y.Z"``
#. ``git tag vX.Y.Z``
#. Test version updated as intended with ``make test-install``


Push to repository
~~~~~~~~~~~~~~~~~~

To do the release, push the tags and the repository state.

#. ``git push``
#. ``git push --tags``

Assuming all the checks pass, this automatically triggers a release on PyPI via the ``.github/workflows/ci-cd-workflow.yml`` action.


Why is there a ``Makefile`` in a pure Python repository?
--------------------------------------------------------

Whilst it may not be standard practice, a ``Makefile`` is a simple way to automate general setup (environment setup in particular).
Hence we have one here which basically acts as a notes file for how to do all those little jobs which we often forget e.g. setting up environments, running tests (and making sure we're in the right environment), building docs, setting up auxillary bits and pieces.

.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _OpenSCM-Runner issue tracker: https://github.com/openscm/openscm-runner/issues
.. _`OpenSCM-Runner's PyPI`: https://pypi.org/project/openscm-runner/
