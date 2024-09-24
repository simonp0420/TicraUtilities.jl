#=
# Introduction
`TicraUtilities` is a package to facilitate working with Ticra-compatible files and data objects in the 
Julia language.
[TICRA](https://www.ticra.com) is a Danish company specializing in antenna analysis and synthesis software. Their
products for analysis and design of reflector antennas (and associated feed structures) are widely regarded as
standards in the antenna community.  As a result, the file formats defined and used by TICRA software have also
become de facto industry standards. 
    
This package provides utilities in the Julia programming language for 
* Reading, writing, and plotting Ticra-compatible cut (.cut) files.
* Reading, parsing, and writing Ticra Object Repository (.tor) files.
* Converting cut files to and from spherical wave expansion (.sph) files.
* Reading and writing so-called "station" (.stn) files.
* Reading and writing array excitation (.exi) files.
* Reading and writing "surface" (.sfc) files.
* And more...
=#
