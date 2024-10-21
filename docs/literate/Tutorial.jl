#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#!nb # # Tutorial
#nb # # Introduction to the TicraUtilities Package for Julia
#nb # [https://github.com/simonp0420/TicraUtilities.jl](https://github.com/simonp0420/TicraUtilities.jl)
#-
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
#nb using Pkg
#nb Pkg.activate(temp=true)
#nb Pkg.add(["TicraUtilities", "Plots"])
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ## Reading and Plotting Cut Files
#
# We begin by reading in an existing cut file:
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
using TicraUtilities

cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "test.cut")
cut = read_cutfile(cutfile)
#!nb # We see that `read_cutfile` returns a `Cut` object.  The information printed to
#!nb # the REPL shows that the `Cut` object contains 72 distinct ϕ cuts, each consisting of 361
#!nb # points in the θ direction, with θ ranging from 0° to 180° in steps of 0.5°. 
#!nb # See [Cut](@ref) for documentation on the remaining fields.
#
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
#nb ?Cut
    
#nb # %% [markdown] {"slideshow": {"slide_type": "fragment"}}
# Values within the `Cut` object can be accessed via the functions
# [`get_ncomp`](@ref), [`get_icut`](@ref), [`get_icomp`](@ref), [`get_text`](@ref), 
# [`get_theta`](@ref), [`get_phi`](@ref), and [`get_evec`](@ref).  Some derived
# quantities can be obtained via the functions [`power`](@ref), [`amplitude_db`](@ref), 
# and [`phase_deg`](@ref).
#
#!nb # This package includes a recipe for plotting a `Cut` object using the standard 
#!nb # [`Plots`](https://docs.juliaplots.org/stable/) package.  So a plot of the cut
#!nb # can be generated very simply:

#nb %% [code] {"slideshow": {"slide_type": "slide"}}
using Plots
plot(cut)

#!nb # The default result in this case is not very satisfactory.  Note that each trace is automatically
#!nb # labeled according to its polarization (either LHCP or RHCP for this case) and ϕ value, but there 
#!nb # are too many ϕ cut labels to fit in the legend.  Let's tweak the plot by just selecting a 
#!nb # just a few ϕ values to plot.  We will also normalize the amplitudes according to the peak 
#!nb # directivity, and print this peak value in the title of the plot.  Finally, we'll adjust the axes
#!nb # limits and tick locations:
#

#nb %% [code] {"slideshow": {"slide_type": "slide"}}
peakdb = round(maximum_db(cut), digits = 3)
plot(
    cut, 
    normalization = :peak,
    phi = 0:90:315, 
    title = "Peak = $peakdb dB",
    framestyle = :box,
    xlim = (0, 180),
    xtick = 0:20:180,
    ylim = (-60, 0),
    ytick = -60:5:0,
)

#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#!nb # The plot appearance is now more pleasing. We used the [`maximum_db`](@ref) function to label the peak
#!nb # field value used for normalizing the plot.  For plotting, most of the keyword arguments that were 
#!nb # passed to the `plot` function above come standard with the `Plots` package, but `phi` and `normalization` 
#!nb # are custom keywords introduced by `TicraUtilities` and are supported
#!nb # only when plotting a `Cut` object.  The full list of custom keywords available 
#!nb # for plotting `Cut` objects is presented in the table below:
#
# ### Custom `Plots` Keywords Specific to `Cut` Objects
# |Keyword        | Description          | Legal values                                  |
# |:--------------|:---------------------|:---------------------------------------------------|
# | phi           | ϕ value(s) to plot   | Scalar or iterable (defaults to [`get_phi(cut)`](@ref)) |
# | theta         | θ value(s) to plot   | Scalar or iterable (defaults to [`get_theta(cut)`](@ref)) |
# | pol           | Polarization to plot | `:both`, `:copol`, `:xpol`, `1`, or `2` (defaults to `:both`) |
# | quantity      | Quantity to plot     | `:db` (the default), `:power`, `:linear`, or `:phase`  |
# | normalization | How to normalize plot | `:peak` or a scalar (defaults to `NaN` meaning don't normalize) |
#
#!nb # Values passed in via the `phi` keyword must all be present in the `Cut` object.
#!nb # If the values passed via the `theta` keyword are spaced more finely than those stored in the `Cut`
#!nb # object, the pattern will be interpolated onto these values using a cubic spline interpolant.  This
#!nb # is demonstrated in the following example plot:

#nb %% [code] {"slideshow": {"slide_type": "slide"}}
scatter(
    cut, 
    label = "Interpolated In θ",
    pol = :copol,
    normalization = :peak,
    phi = 0,
    theta = 0:0.1:5,
    title = "Normalized Copol (LHCP), Peak = $peakdb dB",
    framestyle = :box,
    xlim = (0, 5),
    xtick = 0:0.5:10,
    ylim=(-0.2, 0.01),
    ytick = -10:0.05:0)
scatter!(
    cut, 
    label = "No θ Interpolation",
    pol = :copol,
    normalization = :peak,
    phi = 0,
    marker = :cross,
    ms = 6,
    msw = 2)

#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#
#!nb # The above scatter plot clearly illustrates the results of interpolating in ``\theta``.  The example
#!nb # also shows that the default trace labels appearing in the legend can be overridden as desired.
# ## Manipulation and Conversion of `Cut` Objects
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Symmetric and Asymmetric Cuts
#!nb # In the previous examples the `cut` variable
#!nb # contains "asymmetric" cuts, each beginning at ``\theta = 0``.   A "symmetric" cut would cover equal extents
#!nb # in the positive and negative ``\theta`` directions.  Functions [`asym2sym`](@ref) and [`sym2asym`](@ref) 
#!nb # can be used to convert between these types of cuts.  Continuing with the asymmetric `Cut` object from the
#!nb # previous examples:
scut = asym2sym(cut) # Create a symmetric cut
#
#!nb # Note that the symmetric cuts only extend to 175° in ``\phi``, and that each cut covers the range 
#!nb # ``-180^\circ \le \theta \le 180^\circ``.  We plot the new cut below to see this alternative representation
#!nb # of the same data:
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
plot(
    scut, 
    normalization = :peak,
    phi = 0:45:135,
    title = "Normalized Symmetric Cut, Peak = $peakdb dB",
    framestyle = :box,
    xlim = (-180, 180),
    xtick = -180:30:180,
    ylim = (-60, 0),
    ytick = -60:5:0,
)
#
#!nb # This type of symmetric plot can be useful for spotting pattern asymmetries.
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Changing the Polarization Basis of a `Cut`
# The functions [`convert_cut`](@ref) and [`convert_cut!`](@ref) can be used to change
# which of the following pairs of field components
# 1. ``E_\theta`` and ``E_\phi``
# 2. ``E_R`` and ``E_L`` (CP or Circular Polarization, right- and left-handed, resp.)
# 3. ``E_h`` and ``E_v`` (L3 or Ludwig 3, directed along ``x`` and ``y``, resp.)
# is used for representing the electric field vector stored in a `Cut` object. 
#!nb # The `cut` variable used in the previous examples uses CP components, as can be verified
#!nb # using the [`get_icomp`](@ref) function:

#nb %% [code] {"slideshow": {"slide_type": "slide"}}
get_icomp(cut)
#!nb # The possible values for `icomp` and their meanings are documented in [ttools24; Sec. 9.7, p. 3306](@cite).
#!nb # Alternatively, one can request help on the `Cut` type at the Julia REPL. Note: At present, most of the 
#!nb # functions in `TicraUtilities` support only 1, 2, or 3 as possible values for `icomp`.  In this case, the 
#!nb # value is `2`, meaning that a circular polarization basis was used.
#!nb #
#!nb # If we convert the polarization basis to a Ludwig 3 representation:
#
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cut_L3 = convert_cut(cut, 3)
peakdb_L3 = round(maximum_db(cut_L3), digits = 3)
#!nb # we see that the peak directivity has been reduced by about 3 dB from its previous value.  Since
#!nb # the boresight radiated field is nearly perfectly circularly polarized, then both L3 components
#!nb # will be approximately equal in magnitude as shown in the following plot:
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
plot(
    cut_L3, 
    normalization = :peak,
    phi = 0:90:315, 
    title = "Ludwig 3 Components, Peak = $peakdb_L3 dB",
    framestyle = :box,
    xlim = (0, 180),
    xtick = 0:20:180,
    ylim=(-60, 0),
    ytick = -60:5:0,
)

#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Cut Normalization
#!nb # Typically, the fields recorded in cut files are normalized so that the total radiated power is ``4\pi``.
#!nb # When this is the case, the field magnitude squared is numerically equal to directivity.  The power integral
#!nb # for a `Cut` object can be calculated using the [`power`](@ref) function:

#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
power(cut) / 4π
#!nb # We see that the power in `cut` is very close to ``4\pi``.  Although not really necessary in this case,
#!nb # we can modify `cut` to make the normalization more
#!nb # nearly exact by using the [`normalize!`](@ref) function:
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
normalize!(cut)
power(cut) / 4π
#!nb # After explicitly normalizing the cut, its radiated power is almost exactly equal to ``4\pi``.
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Synthesizing a Cut from E- and H-Plane Patterns
#!nb # A "BOR₁" horn [kildal2015](@cite) is rotationally symmetric and contains only TE₁ₙ and TM₁ₙ waveguide modes in its 
#!nb # radiating aperture.  Its radiated far field can therefore be expressed in terms
#!nb # of the E-plane and H-plane (principal plane) patterns it radiates when excited for linear polarization. 
#!nb # The [`eh2bor1cut`](@ref) function synthesizes a `Cut` object for a BOR₁ antenna from
#!nb # its principal plane patterns, optionally adding in a specified level of crosspol due to an 
#!nb # imperfect feed network.

#!nb # We begin by reading a cut file for a BOR₁ horn created by Ticra's CHAMP program:
#!nb using TicraUtilities
#!nb using Plots
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "ticra_hpol_horn.cut")
cut = read_cutfile(cutfile)
println("Peak =  ", round(maximum_db(cut), digits = 2), " dB")
#
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cut

#-
#!nb # As seen from the above output, the cut file contains 3 cuts at ``\phi = 0^\circ, 45^\circ``, and
#!nb # ``90^\circ``.  Since the cut file was created for horizontal excitation of the horn, the dominant
#!nb # far-field polarization will be Ludwig 3 horizontal, which is stored in the first polarization slot of the cut.
#!nb # It also follows that the E-plane pattern can be extracted from the first (``\phi = 0^\circ``) cut
#!nb # and the H-plane pattern from the last (``\phi = 90^\circ``) cut:
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
fe = get_evec(cut, 1)[:, 1] # ϕ = 0° cut is E-plane for horizontal pol
fh = get_evec(cut, 1)[:, end] # ϕ = 90° cut is H-plane for horizontal pol
plot(xlim=(0, 180),
     xlabel = "θ (deg)",
     ylim = (-60, 30),
     ylabel = "Amplitude (dB)",
     framestyle = :box,
     title = "BOR₁ Horn Principal Plane Patterns",
)
theta = get_theta(cut)
plot!(theta, 10 .* log10.(abs2.(fe)), label = "E-Plane")
plot!(theta, 10 .* log10.(abs2.(fh)), label = "H-Plane")

#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# Suppose now that we wish to create a `Cut` object for this horn, but assuming that it has been
# excited to generate a predominantly RHCP (right-hand circularly polarized) far field:
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cutrhcp = eh2bor1cut(theta, fe, fh; pol = :rhcp)
println("Peak RHCP =  ", round(maximum_db(cutrhcp), digits = 2), " dBi")
println("(Radiated power)/4π = ",  power(cutrhcp) / 4π)
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cutrhcp
#!nb # As shown above, the maximum field magnitude is about the same as that for the principal plane cuts. 
#!nb # And since the cut is properly power-normalized, this value is the maximum partial directivity 
#!nb # to RHCP polarization.
#!nb # As stated in the documentation for [`eh2bor1cut`](@ref), `cutrhcp` contains cuts at 
#!nb # ``\phi = 0^\circ, 90^\circ, 180^\circ, \text{and } 270^\circ``.  Plotting the cut at ``\phi = 0^\circ``:
#-
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
plot(cutrhcp, 
     phi = 0,
     xlim = (0, 90),
     xtick = 0:10:90,
     ylim = (-60, 0),
     framestyle = :box, 
     normalization = :peak,
     title = "BOR₁ Horn Excited for RHCP, Normalized Pattern")
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#!nb # The plot confirms that the dominant polarization is RHCP.  The maximum crosspol level is about 
#!nb # 45 dB below the copol peak.  To simulate the effect of an imperfect feed network that injects 
#!nb # crosspol (LHCP) at a level 30 dB below the desired copol, we can use the `xpd` keyword argument:
#nb # Simulate crosspol injected from imperfect feed network using `xpd`:
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cutrhcp2 = eh2bor1cut(theta, fe, fh; pol = :rhcp, xpd = -30)
plot(cutrhcp2,
     phi = 0,
     xlim = (0, 90),
     xtick = 0:10:90,
     ylim = (-60, 0),
     framestyle = :box,
     normalization = :peak,
     title = "RHCP-Excited BOR₁ Horn, Normalized Pattern\n-30 dB Xpol Added\n")
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "fragment"}}
#!nb # As expected, the boresight crosspol level is now 30 dB below the copol peak, and the crosspol
#!nb # pattern resembles a scaled version of the copol pattern, at least in the vicinity
#!nb # of boresight.  
# Note that the phase of the injected crosspol can be specified using the 
# `xpphase` keyword argument of [`eh2bor1cut`](@ref).
#!nb #
#!nb # After creating the `Cut` object `cutrhcp` in this manner, it can be saved as a Ticra-compatible
#!nb # cut file using [`write_cutfile`](@ref), or it can be converted to a spherical wave expansion using
#!nb # [`cut2sph`](@ref), as discussed below in [SWE/Cut Conversion](@ref).  The resulting `SPHQPartition`
#!nb # object can be saved as a spherical wave expansion (.sph) file, if desired, using [`write_sphfile`](@ref).
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ## Spherical Wave Expansions (SWEs)
# The functions in `TicraUtilities` can read, write, and manipulate Ticra-compatible spherical wave 
# expansion (.sph) files that employ the newer, more accurate, so-called "Q" modes.  To read the 
# contents of such a file, one uses the 
# [`read_sphfile`](@ref) function.  Below we demonstrate the function's use on a spherical wave expension
# file previously generated using GRASP:
#-
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
#!nb using TicraUtilities
sphfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "center_element_rhcp_excited_q.sph")
sph_grasp = read_sphfile(sphfile)
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#!nb # `read_sphfile` returns an object of type [`SPHQPartition`](@ref), containing the SWE data
#!nb # for a single partition (i.e., frequency) stored in a SWE file.  If the file had consisted of 
#!nb # multiple partitions then the returned object would be a vector of such objects.
# A `SPHQPartition` object (or a vector of such objects) can be written to a Ticra-compatible file using
# [`write_sphfile`](@ref).  The values stored in the fields of a `SPHQPartition` object can be retrieved using the 
# functions [`get_prgtag`](@ref), [`get_idstrg`](@ref), [`get_nthe`](@ref), [`get_nphi`](@ref), 
# [`get_nmax`](@ref), [`get_mmax`](@ref), [`get_qsmns`](@ref), and [`get_powerms`](@ref), along with
# [`get_t4`](@ref) through [`get_t8`](@ref).
#
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
#!nb # ### SWE/Cut Conversion
# #### Cut to SWE Conversion
# A `Cut` object, a vector of `Cut` objects, or a cut file can be converted to SWE representation via the [`cut2sph`](@ref)
# function.  
#!nb # We'll demonstrate this conversion using a measured cut file for a central element of a closely spaced
#!nb # array of planar radiating elements.  Because of strong mutual coupling, the pattern is asymmetrical, as shown in
#!nb # the plot below:
#!nb using TicraUtilities
#!nb using Plots
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "center_element_rhcp_excited.cut")
cut_meas = read_cutfile(cutfile)
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
plot(cut_meas, legend=nothing, framestyle=:box)
#!nb # Note that although the cut file contains pattern data out to ``\theta = 180^\circ``, the values for 
#!nb # ``\theta > 90^\circ`` are identically zero.  This fact and the large number
#!nb # of samples in the two angular directions make it an interesting and relatively difficult case for 
#!nb # spherical wave expansion.  We generate the SWE representation as follows:
#-
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
sph_julia = cut2sph(cut_meas)
#
#!nb # On a Core i7-9700 computer running Julia 1.11.0 this conversion takes about 45 msec.  The `sph_grasp` object 
#!nb # read in the previous example was generated by the Ticra GRASP program from the same measured cut file. 
#!nb # The maximum difference between the spherical wave "Q" coefficients in `sph_julia` and `sph_grasp` is
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
maximum(abs, get_qsmns(sph_julia) - get_qsmns(sph_grasp))
# In the next section we show that both expansions provide similar accuracy in reconstructing the far-field pattern.
#-
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# #### SWE to Cut Conversion
#!nb # We can reconstruct the pattern of the previous example by converting the `SPHQPartition` object `sph_julia` 
#!nb # into a new `Cut` object  using the [`sph2cut`](@ref) function:
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
cut_julia = sph2cut(sph_julia) # Round trip via Julia: measured cut -> SWE -> reconstructed cut
#!nb # We will compare this `Cut` object to one from cut file "center\_element\_rhcp\_excited_q.cut".  This latter file
#!nb # is also the result of a round-trip conversion cut ``\rightarrow`` SWE ``\rightarrow`` cut, beginning with the 
#!nb # original measured data cut file, but performed entirely using the TICRA GRASP program. The copol and 
#!nb # crosspol comparison is shown below:
cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "center_element_rhcp_excited_q.cut")
cut_grasp = read_cutfile(cutfile) # Round trip via GRASP
rhcp_err_db =  10 * log10(maximum(abs2, get_evec(cut_julia, 1) - get_evec(cut_grasp, 1)))
lhcp_err_db =  10 * log10(maximum(abs2, get_evec(cut_julia, 2) - get_evec(cut_grasp, 2)))
(rhcp_err_db, lhcp_err_db) # Compare Julia and GRASP cut reconstructions
#-
#!nb # As shown above, the maximum differences between the reconstructed fields via GRASP and Julia functions are 
#!nb # extremely small, on the order of -100 dB.  We can also compare the reconstructed fields versus the original
#!nb # measured cut data:
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
using LinearAlgebra: norm
julia_err_db = 20 * log10(maximum(norm, get_evec(cut_julia) - get_evec(cut_meas)))
grasp_err_db = 20 * log10(maximum(norm, get_evec(cut_grasp) - get_evec(cut_meas)))
(julia_err_db, grasp_err_db) # compare each to orig. measured cut
#-
#!nb # GRASP- and Julia-reconstructed patterns show similar agreement with original measured data.
#!nb # Plots of the reconstructed copol and crosspol patterns confirm the similarity of the GRASP and
#!nb # Julia reconstructed pattern data:
#nb %% [code] {"slideshow": {"slide_type": "slide"}}
#!nb using Plots
plot(xlim = (0, 180), ylim = (-100, 15), framestyle = :box,
     title = "Copol (RHCP) Comparison")
plot!(cut_julia, phi = 0, pol = 1, label = "RHCP sph2cut")
plot!(cut_grasp, phi = 0, pol = 1, ls = :dash, label = "RHCP GRASP")
#-
#nb %% [code] {"slideshow": {"slide_type": "fragment"}}
plot(xlim = (0, 180), ylim = (-100, 15), framestyle = :box,
     title = "Crosspol (LHCP) Comparison")
plot!(cut_julia, phi = 0, pol = 2, label = "LHCP sph2cut")
plot!(cut_grasp, phi = 0, pol = 2, ls = :dash, label = "LHCP GRASP")
#-
#!nb # The rapid dropoff at ``\phi = 90^\circ`` occurs because the original pattern data from which the 
#!nb # spherical wave coefficients were derived was only nonzero in the forward hemisphere of the antenna.
#!nb # Finite amplitudes in the reconstructed patterns for ``\theta > 90^\circ`` are due to the
#!nb # discontinuity in the fields at ``\theta = 90^\circ``, which would require an infinite number of 
#!nb # modes to reproduce exactly.
#!nb # 
#!nb # 
#nb # %% [markdown] {"slideshow": {"slide_type": "slide"}}
# ## Surface Files
# Ticra-compatible surface (.sfc) files can be read using the [`read_surface`](@ref) function, and written to 
# disk using the [`write_surface`](@ref) function.  The results of reading a surface file are stored in an
# object of type [`TicraUtilities.Surface`](@ref) as in the following example:
sfcfile = joinpath(joinpath(dirname(pathof(TicraUtilities)), "..", "test", "parent_parabola.sfc"))
sfc = read_surface(sfcfile)
# Values within the `Surface` object can be accessed via the functions
# [`get_x`](@ref), [`get_y`](@ref), [`get_z`](@ref), and [`get_text`](@ref).
#
# The `+` and `-` operators have been overloaded to work on surfaces, resulting in new surfaces whose
# ``z`` values are the sum or difference of those of the operand surfaces:
sfc2 = sfc + sfc
get_z(sfc2) ≈ 2 * get_z(sfc)
#
sfc3 = sfc - sfc
maximum(abs, get_z(sfc3))

#
# ## Array Excitation Files
# Ticra-compatible array excitation (.exi) files can be read using the [`read_exifile](@ref) function, and written to 
# disk using the [`write_exifile`](@ref) function.  The results of reading an array excitation file are stored in an
# object of type [`Exi`](@ref) as in the following example:
exifile = joinpath(joinpath(dirname(pathof(TicraUtilities)), "..", "test", "beam_A14R.exi"))
exi = read_exifile(exifile)
# Values within the `Exi` object can be accessed via the functions
# [`get_header`](@ref), [`get_ampdb`](@ref) (or [`amplitude_db`](@ref)), 
# [`get_phsdeg`](@ref) (or [`phase_deg`](@ref)), and [`get_ids`](@ref).
# For example:
get_ampdb(exi)
# 
get_phsdeg(exi)
#
# ## Optimization Station Files
# Ticra-compatible optimization station (.sta) files, also known as "Field Directions" file,
# can be read using the [`read_stationfile](@ref) function, 
# and written to disk using the [`write_stationfile`](@ref) function.  The results of reading a station file are 
# stored in a vector of objects of type [`Station`](@ref) as in the following example:
stationfile = joinpath(joinpath(dirname(pathof(TicraUtilities)), "..", "test", "scenario2_coverage.sta"))
stations = read_stationfile(stationfile)
# Values within a `Station` object can be accessed via the functions
# [`get_npoint`](@ref), [`get_u`](@ref), [`get_v`](@ref),  [`get_goal`](@ref),  [`get_weight`](@ref),
#  [`get_ipol`](@ref),  [`get_rot`](@ref),  [`get_att`](@ref),  and [`get_id`](@ref).
#
# ## Ticra Object Repository (TOR) Files
# TOR files can be read and written using the functions [`read_torfile`](@ref) and [`write_torfile`](@ref), 
# respectively.  Here is an example of reading a TOR file:
torfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "tabulated_rim_tor_file.tor")
torobjs = read_torfile(torfile)

# [`read_torfile`](@ref) returns a vector of `TorObj` objects.  Here is the first element of this vector:
torobj = torobjs[1]
# The name and Ticra type of the object are shown, followed by propertynames and their corresponding values.
# These can be extracted from the `TorObj` object using functions [`get_name`](@ref), [`get_objtype`](@ref), 
# [`get_propname`](@ref), and [`get_propval`](@ref). For example:
get_name(torobj)
#
get_objtype(torobj)
#
get_propname(torobj)
#
get_propval(torobj)
# The [`parse_tor_struct`](@ref) function can be used to parse the Ticra struct objects listed in the final
# two elements of `get_propval(torobj)`:
struct1 = parse_tor_struct(get_propval(torobj)[end-1])
#
struct2 = parse_tor_struct(get_propval(torobj)[end])
#
# `struct1` and `struct2` are [`NamedTuple`](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple)s.
# Their field names and values can be obtained using Julia's
# `keys` and `values` functions:
keys(struct2)
#
values(struct2)
# Note that the parsed values of numeric quantities with associated units (such as "x" and "y" in the above example)
# are converted to [`Unitful`](https://github.com/PainterQubits/Unitful.jl) quantities.   For example:
struct2.x
# The purely numeric portion and the units can be extracted using functions supplied by `Unitful`:
using Unitful: unit, ustrip
unit(struct2.x)
#
ustrip(struct2.x)
# Alternatively, both can be converted to strings, if desired:
string(struct2.x)
# Of the units that can occur in a TOR file, "in" is the only one which is modified when translated to Julia.
# In the Julia representation, "inch" is used instead of "in" to avoid confusion with the built-in Julia function 
# [`in`](https://docs.julialang.org/en/v1/base/collections/#Base.in).