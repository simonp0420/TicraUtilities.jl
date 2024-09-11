# # Usage Examples

# ## Reading and Plotting Cut Files
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



# Now plot the resulting `Cut` object, using the default settings:
using Plots
plot(cut)

# We see that