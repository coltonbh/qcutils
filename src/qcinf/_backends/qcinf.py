"""
Cheminformatics algorithms not relying on external libraries.

All functions should expose a qcio interface.
"""

import numpy as np
from qcconst import constants
from qcio import LengthUnit, Structure

from .utils import compute_rmsd, kabsch


def _rmsd(
    struct1: Structure,
    struct2: Structure,
    *,
    align: bool = True,
    units: LengthUnit = LengthUnit.BOHR,
) -> float:
    """Compute the RMSD between two structures.
    Args:
        struct1: The first structure.
        struct2: The second structure.
        align: Whether to align the structures before computing the RMSD. Defaults to
            True.
        units: The units of the RMSD. Defaults to LengthUnit.BOHR.
    Returns:
        The RMSD between the two structures in the specified units.
    """
    rmsd_val = compute_rmsd(
        struct1.geometry,
        struct2.geometry,
        align=align,
    )
    return (
        rmsd_val * constants.BOHR_TO_ANGSTROM
        if units == LengthUnit.ANGSTROM
        else rmsd_val
    )


def _align(struct1: Structure, struct2: Structure) -> Structure:
    """
    Align struct1 onto struct2 using the Kabsch algorithm. Returns a new Structure.

    Args:
        struct1 (Structure): The first structure to align.
        struct2 (Structure): The structure to align to.

    Returns:
        A new Structure object aligned to struct2.
    """
    R, centroid_P, centroid_Q = kabsch(struct1.geometry, struct2.geometry)

    # Rotate and translate struct1 to align with struct2
    aligned_coords1 = np.dot(struct1.geometry - centroid_P, R.T) + centroid_Q

    # Dump the structure to a dictionary and modify the geometry.
    new_struct = struct1.model_dump()
    new_struct["geometry"] = aligned_coords1

    return Structure(**new_struct)
