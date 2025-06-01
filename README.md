# Quantum Chemistry Utilities (qcutils)

[![image](https://img.shields.io/pypi/v/qcutils.svg)](https://pypi.python.org/pypi/qcutils)
[![image](https://img.shields.io/pypi/l/qcutils.svg)](https://pypi.python.org/pypi/qcutils)
[![image](https://img.shields.io/pypi/pyversions/qcutils.svg)](https://pypi.python.org/pypi/qcutils)
[![Actions status](https://github.com/coltonbh/qcutils/workflows/Tests/badge.svg)](https://github.com/coltonbh/qcutils/actions)
[![Actions status](https://github.com/coltonbh/qcutils/workflows/Basic%20Code%20Quality/badge.svg)](https://github.com/coltonbh/qcutils/actions)

Standardized cheminformatics functions such as `rmsd`, structure alignment, `smiles_to_structure`, `structure_to_smiles`, and many more using [qcio](https://qcio.coltonhicks.com/) data structures.

`qcutils` works in harmony with a suite of other quantum chemistry tools for fast, structured, and interoperable quantum chemistry.

## The QC Suite of Programs

- [qcconst](https://github.com/coltonbh/qcconst) - Physical constants, conversion factors, and a periodic table with clear source information for every value.
- [qcio](https://github.com/coltonbh/qcio) - Elegant and intuitive data structures for quantum chemistry, featuring seamless Jupyter Notebook visualizations. [Documentation](https://qcio.coltonhicks.com)
  [qcutils](https://github.com/coltonbh/qcutils) - Standardized cheminformatics functions such as `rmsd`, structure alignment, `smiles_to_structure`, `structure_to_smiles`, and many more using [qcio](https://qcio.coltonhicks.com/) data structures.
- [qcparse](https://github.com/coltonbh/qcparse) - A library for efficient parsing of quantum chemistry data into structured `qcio` objects.
- [qcop](https://github.com/coltonbh/qcop) - A package for operating quantum chemistry programs using `qcio` standardized data structures. Compatible with `TeraChem`, `psi4`, `QChem`, `NWChem`, `ORCA`, `Molpro`, `geomeTRIC` and many more.
- [BigChem](https://github.com/mtzgroup/bigchem) - A distributed application for running quantum chemistry calculations at scale across clusters of computers or the cloud. Bring multi-node scaling to your favorite quantum chemistry program.
- `ChemCloud` - A [web application](https://github.com/mtzgroup/chemcloud-server) and associated [Python client](https://github.com/mtzgroup/chemcloud-client) for exposing a BigChem cluster securely over the internet.

## Installation

```bash
python -m pip install qcutils
```

## Support

If you have any issues with `qcutils` or would like to request a feature, please open an [issue](https://github.com/coltonbh/qcutils/issues).
