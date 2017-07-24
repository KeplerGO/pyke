==========
PyKE tasks
==========

Calling tasks from Python
-------------------------

PyKE provides a set of tasks which can be called both as Python functions
and as command-line utilities. For example, PyKE provides the ``kepextract`` task which enables a lightcurve to be extracted from a Target Pixel File.  You can import and execute the task in Python using::

    >>> from pyke import kepextract
    >>> kepextract("target-pixel-file.fits", "output-lightcurve.fits")

If you are using an interactive IPython shell, you can learn about the use of the task and its optional parameters by adding a question mark behind the function name and hitting return::

    >>> kepextract?


Calling tasks from a terminal
-----------------------------

For convenience, many of PyKE's tasks are installed as command-line utilities on your system during installation.  This allows you to perform common tasks from your favorite terminal without having to open a Python interpreter. For example, the above task can be executed using::

    kepextract target-pixel-file.fits output-lightcurve.fits 

And its help documentation can can be retrieved using::

    kepextract --help

More generally, most of PyKE's tasks can be run from the command-line using the following pattern::

    name-of-the-tool input-file output-file --options

where ``input-file`` is either a light curve or a target pixel file and
``output-file`` is the destination of the output.  Click on a task in the list below to retrieve detailed documentation.


List of all PyKE tasks
----------------------

.. toctree::
    :maxdepth: 2

    api/index
