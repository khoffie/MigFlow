using Pkg
Pkg.activate(".")
using StatsBase, Serialization, Turing, Plots
germ = loadallGermData(; sample = false)
us = loadallUSdata()

savelp = function()
    path = "./manuscript_input"
    out = Dict{String, Any}()  # Create a dictionary to store the plots
    ages = ["below18", "18-25", "25-30", "30-50", "50-65", "above65"]
    for age in ages
        chain = deserialize(path * "/germchain_$(age)")
        out[age] = Plots.plot(chain[:lp], title = "LP for $(age)")
    end
    p = Plots.plot(out["below18"], out["18-25"], out["25-30"],
                   out["30-50"], out["50-65"], out["above65"],
                   layout = (3, 2))
    savefig(path * "/plots/lp_values.pdf")
end


savelp()


