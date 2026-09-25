"""ddatoolbox: shared library for the ShinMassBank DDA annotation pipeline.

Bundles the sub-modules used by the top-level ``DDA*.py`` pipeline scripts:

- ``ddatoolbox``: spectral clique construction, averaging, and MS1/hive file I/O.
- ``massbankio``: the ``MassBankRecord`` class and MassBank flat-file record I/O.
- ``massbankserver``: local MassBank/GNPS spectral-library search server.
- ``mfsearch``: molecular-formula search (Seven Golden Rules) and mass calculators.
- ``compoundsearch``: wrappers for external annotation tools (SIRIUS, MS-FINDER,
  CFM-ID, MetFrag) and spectral-similarity search.
- ``paths``: dataset/resource folder-path resolution and settings loading.
- ``resources``: cached lookup tables (InChIKey-to-ClassyFire, PubChem-checked list).
"""
from . import ddatoolbox
from . import massbankio
from . import massbankserver
from . import mfsearch
from . import compoundsearch
from . import paths
from . import resources






