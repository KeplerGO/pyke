============
Contributing
============

PyKE is actively developed on its `GitHub repository <https://github.com/KeplerGO/PyKE>`_.

If you encounter a problem with PyKE, we encourage you to `open an issue on the GitHub repository <https://github.com/KeplerGO/PyKE/issues>`_
or to e-mail the Kepler/K2 office at keplergo@mail.arc.nasa.gov.
We also welcome suggestions for enhancements or new features in this way.

If you would like to contribute a patch for a bugfix, please go ahead and open a pull request.

If you want to make a more significant change, such as adding a new tool or feature,
we recommend discussing the change in advance by opening a GitHub issue.


Proposing a change or bugfix to PyKE using a GitHub pull request
----------------------------------------------------------------

If in doubt on how to open a pull request, we recommend Astropy's
"`How to make a code contribution <http://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_" tutorial.
An example is shown below.

1. Fork KeplerGO/PyKE by going to https://github.com/KeplerGO/PyKE and clicking on "Fork"

2. Clone your fork to your computer:
    $ git clone https://github.com/your-gh-handle/PyKE.git

3. Install the development version of PyKE:
    $ cd PyKE
    $ pip install -e .

4. Add the KeplerGO remote to your github enviroment:
    $ git remote add upstream https://github.com/KeplerGO/PyKE.git

5. Let's make sure everything is good, do:
    $ git remote -v

You should see something like:
    origin	https://github.com/mirca/PyKE.git (fetch)
    origin	https://github.com/mirca/PyKE.git (push)
    upstream	https://github.com/KeplerGO/PyKE.git (fetch)
    upstream	https://github.com/KeplerGO/PyKE.git (push)

6. Now you are ready to start contributing; make a branch and check out to it
    $ git branch name-of-my-branch
    $ git checkout name-of-my-branch

7. Do the changes you want and add them:
    $ git add .

8. Commit and push your changes:
    $ git commit -m "that's my changes y'all"
    $ git push origin name-of-my-branch

9. That's it!! =) Go to https://github.com/KeplerGO/PyKE and you should see a button "Compare and open a pull request" =)
