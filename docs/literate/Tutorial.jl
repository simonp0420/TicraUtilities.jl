# # Tutorial

# ## Reading and Plotting Cut Files
#
# We begin by reading in an existing cut file:
using TicraUtilities

cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "test.cut")
cut = read_cutfile(cutfile)
#
# We see that `read_cutfile` returns a `Cut` object.  The information printed to
# the REPL shows that the `Cut` object contains 72 distinct ϕ cuts, each consisting of 361
# points in the θ direction, with θ ranging from 0° to 180° in steps of 0.5°. 
# See [Cut](@ref) for documentation on the remaining fields.
#
# Values within the `Cut` object can be accessed via the functions
# [`get_ncomp`](@ref), [`get_icut`](@ref), [`get_icomp`](@ref), [`get_text`](@ref), 
# [`get_theta`](@ref), [`get_phi`](@ref), and [`get_evec`](@ref).  Some derived
# quantities can be obtained via the functions [`power`](@ref), [`amplitude_db`](@ref), 
# and [`phase_deg`](@ref).
#
# This package includes a recipe for plotting a `Cut` object using the standard 
# [`Plots`](https://docs.juliaplots.org/stable/) package.  So a plot of the cut
# can be generated very simply:
using Plots
plot(cut)

# The default result in this case is not very satisfactory.  Note that each trace is automatically
# labeled according to its polarization (either LHCP or RHCP for this case) and ϕ value, but there 
# are too many ϕ cut labels to fit in the legend.  Let's tweak the plot by just selecting a 
# just a few ϕ values to plot.  We will also normalize the amplitudes according to the peak 
# directivity, and print this peak value in the title of the plot.  Finally, we'll adjust the axes
# limits and tick locations:

peakdb = round(maximum_db(cut), digits=3)
plot(
    cut, 
    normalization = :peak,
    phi=0:90:315, 
    title="Peak = $peakdb dB",
    framestyle=:box,
    xtick=0:20:180, xlim=(0,180),
    ytick=-60:5:0, ylim=(-60,0),
)

# The plot appearance is now more pleasing. We used the [`maximum_db`](@ref) function to label the peak
# field value used for normalizing the plot.  For plotting, most of the keyword arguments that were 
# passed to the `plot` function above come standard with the `Plots` package, but `phi` and `normalization` 
# are custom keywords introduced by `TicraUtilities` and are supported
# only when plotting a `Cut` object.  The full list of custom keywords available 
# for plotting `Cut` objects is presented in the table below:
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
# Values passed in via the `phi` keyword must all be present in the `Cut` object.
# If the values passed via the `theta` keyword are spaced more finely than those stored in the `Cut`
# object, the pattern will be interpolated onto these values using a cubic spline interpolant.  This
# is demonstrated in the following example plot:

scatter(
    cut, 
    label = "Interpolated In θ",
    pol = :copol,
    normalization = :peak,
    phi=0,
    theta=0:0.1:5,
    title="Normalized Copol (LHCP), Peak = $peakdb dB",
    framestyle=:box,
    xtick=0:0.5:10, xlim=(0,5),
    ytick=-10:0.05:0, ylim=(-0.2,0.01),
)
scatter!(
    cut, 
    label = "No θ Interpolation",
    pol = :copol,
    normalization = :peak,
    phi=0,
    marker = :cross, ms = 6, msw = 2,
)

# The above scatter plot clearly illustrates the results of interpolating in ``\theta``.  The example
# also shows that the default trace labels appearing in the legend can be overridden as desired.
#
# ## Manipulation and Conversion of `Cut` Objects
# ### Symmetrical and Asymmetrical Cuts
# In the previous examples the `cut` variable
# contains "asymmetric" cuts, each beginning at ``\theta = 0``.   A "symmetric" cut would cover equal extents
# in the positive and negative ``\theta`` directions.  Functions [`asym2sym`](@ref) and [`sym2asym`](@ref) 
# can be used to convert between these types of cuts.  Continuing with the asymmetric `Cut` object from the
# previous examples:
scut = asym2sym(cut) # Create a symmetric cut
#
# Note that the symmetric cuts only extend to 175° in ``\phi``, and that each cut covers the range 
# ``-180^\circ \le \theta \le 180^\circ``.  We plot the new cut below to see this alternative representation
# of the same data:
plot(
    scut, 
    normalization = :peak,
    phi=0:45:135,
    title="Normalized Symmetric Cut, Peak = $peakdb dB",
    framestyle=:box,
    xtick=-180:30:180, xlim=(-180,180),
    ytick=-60:5:0, ylim=(-60,0),
)
# This type of symmetric plot can be useful for spotting pattern asymmetries.
#
# ### Changing the Polarization Basis of a `Cut`
# The functions [`convert_cut`](@ref) and [`convert_cut!`](@ref) can be used to change
# which of the following pairs of field components
# 1. ``E_\theta`` and ``E_\phi``
# 2. ``E_R`` and ``E_L`` (CP or Circular Polarization, right- and left-handed, resp.)
# 3. ``E_h`` and ``E_v`` (L3 or Ludwig 3, directed along ``x`` and ``y``, resp.)
# is used for representing the electric field vector stored in a `Cut` object. The `cut`
# variable used in the previous examples uses CP components, as can be verified
# using the [`get_icomp`](@ref) function:
get_icomp(cut)
#
# If we convert to an L3 representation:
cut_L3 = convert_cut(cut, 3)
peakdb_L3 = round(maximum_db(cut_L3), digits=3)
# we see that the peak directivity has been reduced by about 3 dB from its previous value.  Since
# the boresight radiated field is nearly perfectly circularly polarized, then both L3 components
# will be approximately equal in magnitude as shown in the following plot:
plot(
    cut_L3, 
    normalization = :peak,
    phi=0:90:315, 
    title="Ludwig 3 Components, Peak = $peakdb_L3 dB",
    framestyle=:box,
    xtick=0:20:180, xlim=(0,180),
    ytick=-60:5:0, ylim=(-60,0),
)

# ### Cut Normalization
# Typically, the fields recorded in cut files are normalized so that the total radiated power is ``4\pi``.
# When this is the case, the field magnitude squared is numerically equal to directivity.  The power integral
# for a `Cut` object can be calculated using the `power` function:
power(cut) / 4π
# We see that the power in `cut` is very close to ``4\pi``.  We can modify `cut` to make the normalization more
# nearly exact by using the `normalize!` function:
normalize!(cut)
power(cut) / 4π
# The remaining departure of the radiated power from exact equality to ``4\pi`` is due to floating point error.

# ### Synthesizing a Cut from E- and H-Plane Patterns
# A "BOR₁" horn [kildal2015](@cite) is rotationally symmetric and contains only TE₁ₙ and TM₁ₙ waveguide modes in its 
# radiating aperture.  It's radiated far field can therefore be expressed in terms
# of the E-plane and H-plane (principal plane) patterns it radiates when excited for linear polarization. 
# The [`eh2bor1cut`](@ref) function can be used to synthesize a `Cut` object from its principal plane
# patterns, optionally adding in a specified level of crosspol due to an imperfect feed network.

# We begin by reading a cut file for a BOR₁ horn created by Ticra's CHAMP program:
using TicraUtilities
using Plots
cutfile = joinpath(dirname(pathof(TicraUtilities)), "..", "test", "ticra_hpol_horn.cut")
cut = read_cutfile(cutfile)

# We see from the above output that the cut file contains 3 cuts at ``\phi = 0^\circ, 45^\circ``, and
# ``90^\circ``.  Since the cut file was created for horizontal excitation of the horn, the dominant
# far-field polarization will be Ludwig 3 horizontal, which is stored in the first polarization slot of the cut.
# It also follows that the E-plane pattern can be extracted from the first (``\phi = 0^\circ``) cut
# and the H-plane pattern from the last (``\phi = 90^\circ``) cut:
fe = get_evec(cut, 1)[:, 1] # ϕ = 0° cut is E-plane for horizontal pol
fh = get_evec(cut, 1)[:, end] # ϕ = 90° cut is H-plane for horizontal pol
plot(xlim=(0,180),
     ylim=(-60,30),
     framestyle=:box,
     title = "BOR₁ Horn Principal Plane Patterns",
     xlabel="θ (deg)",
     ylabel="Amplitude (dB)")
theta = get_theta(cut)
plot!(theta, 10 .* log10.(abs2.(fe)), label="E-Plane")
plot!(theta, 10 .* log10.(abs2.(fh)), label="H-Plane")

# Suppose now that we wish to create a `Cut` object for this horn, but assuming that it has been
# excited to generate a predominantly RHCP (right-hand circularly polarized) far field.
cutrhcp = eh2bor1cut(theta, fe, fh; pol=:rhcp)
# `cutrhcp` contains cuts at ``\phi = 0^\circ, 90^\circ, 180^\circ, \text{and } 270^\circ``.

plot(cutrhcp, phi=0, xlim=(0,90), ylim=(-50,0), framestyle=:box, normalization=:peak)

# The plot confirms that the dominant polarization is RHCP.  The maximum crosspol level is about 
# 45 dB below the copol peak.  To simulate the effect of an imperfect feed network that injects 
# crosspol (LHCP) at 30 dB below the copol level, we can use the `xpd` keyword argument:
cutrhcp2 = eh2bor1cut(theta, fe, fh; pol=:rhcp, xpd=-30)
plot(cutrhcp2, phi=0, xlim=(0,90), ylim=(-50,0), framestyle=:box, normalization=:peak)
# As expected, the boresight crosspol level is now 30 dB below the copol peak, and the crosspol
# pattern resembles a scaled version of the copol pattern, at least in the vicinity
# of boresight.  Note that the phase of the injected crosspol can be specified using the 
# `xpphase` keyword argument of [`eh2bor1cut`](@ref).

# ## Spherical Wave Expansions