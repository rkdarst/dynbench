Dynamic benchmark resources
===========================

This code and documentation supports the paper "Benchmark model to
assess community structure in evolving networks" [1]_.

For full information on the code, please see ``MANUAL.tex``.  A PDF
version is available under `releases <https://github.com/rkdarst/dynbench/releases>`_.

Dependencies:

* `networkx <http://networkx.github.io/>`_: standard package for
  network analysis.  Found in most Linux distribution repositories.
  Installation is standard with ``setup.py`` or ``pip``.
* `pcd <https://git.becs.aalto.fi/rkdarst/pcd>`_: My own add-on
  package for network and community analysis.  Not in any
  distribution, but can simply be placed in ``PYTHONPATH``.

For convenience, I have a ``.tar.bz2`` that has the benchmark code,
plus local copies of ``networkx`` and ``pcd`` included.  This is the
``dynbench-fulldist-*`` from the corresponding `github releases
<https://github.com/rkdarst/dynbench/releases>`_ page.



Version history
---------------
* The original version corresponding to our first paper [1]_ is at tag
  ``v0.2.0``.  Downloads are at
  https://github.com/rkdarst/dynbench/releases/v0.2.0.  At this time,
  future versions maintain compatibility and default options, so the
  latest version can be used instead.

References
----------

.. [1] Clara Granell, Richard K. Darst, Alex Arenas, Santo Fortunato, and
   Sergio Gómez.  *Benchmark model to assess community structure in
   evolving networks* Phys. Rev. E **92**, 012805.  Published 10 July 2015

   Links: `PRE <http://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012805>`_, `arxiv <http://arxiv.org/abs/1501.05808>`_
