.. _installation:

************
Installation
************

Requirements
============

PyKE has the following requirements, all of which tend to be available by default in a modern installation of Python:

- Python: 2.7, 3.5, 3.6 or later.
- Astropy: 1.3 or later.
- Numpy: 1.11.3 or later.
- Scipy: 0.17 or later.
- Matplotlib: 1.5.3 or later.

If you do not have Python installed, we recommend using the `Anaconda Python <https://www.continuum.io/downloads>`_ distribution, which will install Python in your user directory alongside its most common scientific packages, including all those listed above.



Installing PyKE
===============

Using pip
---------

The easiest way to install or upgrade PyKE is with ``pip``, simply run::

    pip install pyketools --upgrade --no-deps


.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

.. note::

    If you get a ``PermissionError`` this means that you do not have the
    required administrative access to install new packages to your Python
    installation.  In this case you may consider using the ``--user`` option
    to install the package into your home directory.  You can read more
    about how to do this in the `pip documentation
    <http://www.pip-installer.org/en/1.2.1/other-tools.html#using-pip-with-the-user-scheme>`_.


Using the development version
-----------------------------

Alternatively, if you want to experiment with the latest development version of
PyKE, you can install it straight from GitHub::

    $ git clone https://github.com/KeplerGO/PyKE.git
    $ cd PyKE
    $ pip install -e .


Building documentation
======================

.. note::

    Building the documentation is in general not necessary unless you
    are writing new documentation or do not have internet access, because
    the latest (and archive) versions of PyKE's documentation should
    be available at `pyke.readthedocs.io <http://pyke.readthedocs.io>`_ .

To build the documentation, you can do::

    cd docs
    make html

The documentation will be built in the ``docs/_build/html`` directory, and can
be read by pointing a web browser to ``docs/_build/html/index.html``.
