# TicraUtilities - Tools for working with TICRA-compatible files and data objects in the Julia language


| **Documentation**   |  **Tests**     | **CodeCov**  |
|:--------:|:---------------:|:-------:|
|[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simonp0420.github.io/TicraUtilities.jl/stable)  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simonp0420.github.io/TicraUtilities.jl/dev) | [![Build Status](https://github.com/simonp0420/TicraUtilities.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simonp0420/TicraUtilities.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![codecov.io](https://codecov.io/github/simonp0420/TicraUtilities.jl/coverage.svg?branch=main)](https://codecov.io/github/simonp0420/TicraUtilities.jl?branch=main) |




`TicraUtilities` is a package to facilitate working with TICRA-compatible files and data objects in the Julia language. [TICRA](https://www.ticra.com) is a Danish company specializing in antenna analysis and synthesis software. Their products for analysis and design of reflector antennas (and associated feed structures) are widely regarded as standards in the antenna community. As a result, the file formats defined and used by TICRA software have also become de facto industry standards.

This package provides utilities in the Julia programming language for
* Reading, writing, and plotting TICRA-compatible cut (.cut) files.
* Reading and writing spherical wave expansion (.sph) files.
* Converting cut files to and from spherical wave expansion files.
* Reading and writing tabulated electrical properties (TEP) files, both
  the original scattering surface format introduced in GRASP8, and the newer
  periodic unit cell format, introduced in QUPES.
* Bidirectional conversion between old and new TEP types.
* Reading, parsing, and writing TICRA Object Repository (.tor) files.
* Reading and writing so-called "station" (.stn) files.
* Reading and writing array excitation (.exi) files.
* Reading and writing "surface" (.sfc) files.
* And more...



## Package Installation
You can obtain TicraUtilities using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):

```julia
(v1.10) pkg> add TicraUtilities
```

(and then hitting `<Backspace>` to return to the REPL), or with `using Pkg; Pkg.add("TicraUtilities")`.

## Documentation
- Check out the [tutorial video](https://youtu.be/rhinZ-MO9Q4?si=k9A_YCdcqahB5c0J)
- The user manual for the current release is [here](https://simonp0420.github.io/TicraUtilities.jl/stable)
- The user manual for the the development version is [here](https://simonp0420.github.io/TicraUtilities.jl/dev)

## Community
* If there are features you'd like to see added, or if you have other suggestions or questions, please open an 
  [issue](https://github.com/simonp0420/TicraUtilities.jl/issues).
* Pull Requests (PRs) are also welcome!

## News
* Version 1.1: Added functions to read and write TEP files, and convert between old and new TEP file formats.
* Version 1.2: Added functions `issym` and `symsqueeze`.  Several bug fixes.
* Version 1.3: Added function `discard_pol3`. `cut2sph` now accepts symmetric cuts and cuts with 3 polarization slots.

