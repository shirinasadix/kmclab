# kmclab

A lightweight Python package for lattice-based kinetic Monte Carlo (KMC) utilities.  
Currently includes simple hexagonal and square lattice KMC classes.

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/shirinasadix/kmclab.git
cd kmclab
pip install -e .
```

## Usage

```from kmclab import hex_kmc, square_kmc

hex_model = hex_kmc(...)
sq_model = square_kmc(...)
```

## Package Structure 

kmclab/
├── src/
│   └── kmclab/
│       ├── hex.py
│       ├── square.py
│       └── __init__.py
├── tests/
├── pyproject.toml
└── README.md

## License

```
XXX
```
