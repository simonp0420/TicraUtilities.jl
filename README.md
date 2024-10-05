# TicraUtilities
[TICRA](https://www.ticra.com) is a Danish company specializing in antenna analysis and synthesis software. Their
products for analysis and design of reflector antennas (and associated feed structures) are widely regarded as
standards in the aerospace industry.  As a result, the file formats defined and used by TICRA software have also
become de facto industry standards.

This package provides utilities in the Julia programming language for
* Reading, writing, and plotting Ticra-compatible cut (.cut) files.
* Reading, parsing, and writing Ticra Object Repository (.tor) files.
* Converting cut files to and from spherical wave expansion (.sph) files.
* Reading and writing so-called "station" (.stn) files.
* Reading and writing array excitation (.exi) files.
* Reading and writing "surface" (.sfc) files.
* And more...


This package is still under development.  To do:
- [x] Additional functionality
- [x] Improve docstrings
- [ ] Formal documentation

 
[![Build Status](https://github.com/simonp0420/TicraUtilities.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simonp0420/TicraUtilities.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Package Installation
You can obtain TicraUtilities using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):

```julia
(v1.10) pkg> add TicraUtilities
```

or with `using Pkg; Pkg.add("TicraUtilities")`.
