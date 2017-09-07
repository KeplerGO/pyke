.. _tasks:

==========
PyKE tasks
==========

Calling tasks from Python
-------------------------

PyKE provides a set of tasks which can be called both as Python functions
and as command-line utilities. For example, PyKE provides ``kepextract``,
which creates a lightcurve from a Target Pixel File.
You can import and use ``kepextract`` in Python as follows::

    >>> from pyke import kepextract
    >>> kepextract("target-pixel-file.fits", "output-lightcurve.fits")

This will create a file called ``output-lightcurve.fits`` on your current
directory, which stores the extracted lightcurve. Subsequently, this lightcurve
can be visualized using ``kepdraw``::

    >>> from pyke import kepdraw
    >>> kepdraw("output-lightcurve.fits")

If you are using an interactive IPython shell, you can learn about the use
of every task and its optional parameters by adding a question mark behind
the function name and hitting return::

    In [1]: kepextract?


Calling tasks from a terminal
-----------------------------

For convenience, many of PyKE's tasks are installed as command-line utilities
on your system. This allows you to perform common tasks from your favorite
terminal without having to open a Python interpreter. For instance, the example
showed in the previous section can be done from a terminal session as follows::

    $ kepextract target-pixel-file.fits output-lightcurve.fits

And its help documentation can can be retrieved using::

    $ kepextract --help

More generally, most of PyKE's tasks can be run from a terminal using the following pattern::

    $ name-of-the-tool input-file output-file --options

where ``input-file`` is either a light curve or a target pixel file and
``output-file`` is the destination of the output.
Click on a task in the list below to retrieve detailed documentation.


.. toctree::
    :maxdepth: 2

    api/index
