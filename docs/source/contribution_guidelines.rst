How to contribute
=================

We welcome contributions!

You can help by:

* Suggesting new features
* Reporting/fixing bugs
* Adding features to the codebase
* Expanding the testing suit
* Improving the documentation



How to Report Issues
---------------------

When reporting issues, please include as many of the the following details as possible:

* Which version of Kadlu you are using
* The source code that generated the problem (if applicable)
* Which platform you are using (Operating system and version)
* A minimum example that reproduces the issue
* What result you got
* What you were expecting

We use GitLab as a repository, so the best way to submit issues is through the `issues system <https://gitlab.meridian.cs.dal.ca/public_projects/kadlu/issues>`_.

Workflow for Merge Requests
----------------------------

In order make a contribution to the repository, please fork from the ``master`` branch and make a clone of your fork on your local machine.
Make your changes locally. When you are done, commit them to your fork and create a merge request detailing what you did, as well as the reasons why.

Similarly, your commit messages should briefly detail the changes you have made, for example:

.. code-block:: bash

    git commit example.py -m "added example to the docstring of the Ocean::bathy method"


If you are writing a new feature, please ensure you write appropriate test cases and place them under ``kadlu/tests/``.
There are numerous fixtures you can use in ``conftest.py`` and tes/assets contains files that can be used for tests. It is better to use what is already there before adding new fixtures and/or files.

If yours tests need to create temporary files, place them under "tests/assets/tmp". This directory is cleaned by our continous integration setup after the tests are run.

Finally, please run *all* the tests and ensure that *all* the tests complete locally before submitting a merge request.



Thank you for your help!


Running the tests
-----------------

*Kadlu* includes a battery of tests. They are included in the /kadlu/tests/  directory.
We use pytest and doctests.

To run all tests go to the base of your directory

.. code-block:: bash

    cd kadlu-clone
    ls
    
    docker  docs  kadlu  LICENSE.txt  README.md  config.ini  environment.yml  install_dep.sh  requirements.txt  setup.py


and run: ::

    pytest --doctest-modules

You can also specify a module: ::

    pytest kadlu/tests/geospatial/test_interpolation.py
