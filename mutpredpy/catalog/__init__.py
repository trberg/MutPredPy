"""
This module imports the Catalog class and makes it available for external use.

It defines the `__all__` variable to explicitly declare that only the `Catalog` 
class should be exposed when using `from module import *`.
"""

from .catalog import Catalog, CatalogJob

__all__ = ["Catalog", "CatalogJob"]
