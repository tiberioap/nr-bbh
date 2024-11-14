# **Numerical Relativistic Binary Black Hole**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14166347.svg)](https://doi.org/10.5281/zenodo.14166347)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/tiberioap/grav_waldo/blob/main/LICENSE)
[![Orcid](https://img.shields.io/badge/orcid-A6CE39?style=flat&logo=orcid&logoColor=white)](https://orcid.org/0000-0003-1856-6881)

This project implements a numerical simulation of binary black hole (BBH) collisions based on my master's thesis, "[3+1 Formalism of General Relativity](https://repositorio.ufrn.br/handle/123456789/25308)." Written in C, this code simulates BBH dynamics, utilizing Gnuplot to visualize scalar, vector, and tensor fields.

## How It Works
The code leverages the **3+1 formalism** of general relativity, specifically using the BSSN formulation with customized initial data. Key features include:
- **Gauge Conditions**:
  - [1+log](http://doi.org/10.1103/physrevd.52.2059) slicing for the lapse function ($\alpha$);
  - [Gamma-driver](http://doi.org/10.1103/physrevd.73.124011) shift condition for the shift vector ($\beta^i$);
- **BSSN Formulation**: Follows the [Shibata-Nakamura](http://doi.org/10.1103/physrevd.52.5428) and [Baumgarte-Shapiro](http://doi.org/10.1103/physrevd.59.024007) (BSSN) approach;
- **Conformal Factor**: Defined as $\chi = \psi^{-4} = e^{-4\phi}$ ([Campanelli et al.](http://doi.org/10.1103/physrevlett.96.111101));
- **Initial Data**:
  - [Brill-Lindquist](http://doi.org./10.1063/1.1704020) data for **head-on collisions**;
  - [Bowen-York](http://doi.org/10.1103/physrevd.21.2047) data for **quasi-circular orbits**.

### Numerical Methods
- **Evolution Equations**: Solved by [Iterative Crank-Nicholson](http://doi.org/10.1103/physrevd.61.087501) method;
- **Stabilization**: Uses [Kreiss-Oliger](https://www.semanticscholar.org/paper/Methods-for-the-approximate-solution-of-time-Kreiss/283319b5fd1f578b66d40db8e26fa9f587bbd396) dissipation for numerical stability.

All initial conditions and simulation parameters can be modified in `params.h`. Run `run.sh` to execute the simulation and generate visualizations.

For further mathematical and implementation details, refer to `notes.pdf`.

## Limitations
This implementation does **not** include:
- Adaptive Mesh Refinement (AMR);
- Mass and angular momentum extraction;
- Gravitational wave extraction.

**Note**: The AMR method dynamically refines the grid based on spatial curvature, which is critical for enhanced accuracy in extracting mass, angular momentum, and gravitational waves in BBH systems.

## Requirements
To compile and run this code, ensure the following dependencies are installed:
- [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library);
- [FFTW3](http://www.fftw.org/) (Fastest Fourier Transform in the West);
- [Gnuplot](http://www.gnuplot.info/) for visualization.

## BibTeX Citation

```
@software{14166347,
  author       = {Pereira, Tibério},
  title        = {tiberioap/nr-bbh: First Release},
  month        = nov,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.14166347},
  url          = {https://doi.org/10.5281/zenodo.14166347}
}
```