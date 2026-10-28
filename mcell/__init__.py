"""mcell -- analysis code of the multicellularity comparative-genomics project.

Modules
-------
hogs      Streaming OrthoXML (Ancestra HOGs) -> ancestral gene sets and per-branch
          gain / loss / duplication tables.
simtruth  Independent per-branch event counts from the true gene trees of
          ``ancestra.simulate`` (used to validate :mod:`mcell.hogs`).
"""

__version__ = "0.1.0"
