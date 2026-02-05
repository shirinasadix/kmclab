# kmclab

A lightweight Python package for lattice-based kinetic Monte Carlo (KMC) utilities.  
Currently includes simple hexagonal and square lattice KMC classes.

## Installation

Install using pip
```
pip install kmclab
```

To upgrade to the latest version
```
pip install -U kmclab
```

Or you may clone the repository and install in editable mode:

```bash
git clone https://github.com/shirinasadix/kmclab.git
cd kmclab
pip install -e .
```

## Usage

```
from kmclab import hex_kmc

GoombKMC = hex_kmc(n_atoms = 5, n_defects = 5, n_adsorbates = 2, lattice_size = 10, n_steps = 50)

GoombKMC.run()

GoombKMC.anim1panels(filename = 'wtf1')

GoombKMC.anim2panels(filename = 'ov10hy0')

GoombKMC.msdplot(filename = 'MSD_Trajectory')
```

## Package Structure 

```
kmclab/
├── src/
│   └── kmclab/
│       ├── hex.py
│       ├── square.py
│       └── __init__.py
├── tests/
├── pyproject.toml
└── README.md
```

## Demo 

![demo](wtf1.gif)


## License

```
XXX
```




