# O2Rivet project

This directory contains the O2Rivet project.  The aim of the project
is to be able to read in existing or future simulation data and
process it with Rivet analyses.

This has several other benefits

- We can leverage already existing simulation productions.
- Using standard ALICE grid tools, users can quickly and easily run
  their Rivet analyses on large data samples.
- Event generators that are available to the O2 simulations are also
  available to Rivet.
- Development of Rivet analyses becomes easier since they integrate
  into the regular O2 pipe-line.
- Robust comparison of results to models predictions, etc., becomes
  easier.
- Event generators can be validated against known results as part of
  the production.

The Data-Processing-Layer (DPL) `o2-analysis-mm-rivet` can be run as
separate jobs or as part of a larger pipeline.

## Content of this package

- [`doc`](doc): Further documentation and examples.
- [`DataModel`](DataModel): Output data structure
  `o2::rivet::RivetAOs`
- [`Tasks`](Tasks): The DPL `o2-analysis-mm-rivet`
- [`macros`](macros):  Some utility ROOT scripts
- [`Tools`](Tools): Other utility scripts

<!-- EOF -->

