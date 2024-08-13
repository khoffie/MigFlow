gen_preds = function(mapmodel, optis)
    ## expects values in col2, names in col1
    chain = Chains([optis[: , 2]], optis[: , 1]) 
    preds = generated_quantities(mapmodel, chain)
    return preds
end

get_params = function(mapfit) 
    params = DataFrame(names = names(mapfit.values, 1), values = mapfit.values.array)
    return params   
end

load_flows = function()
    if ENV["USER"] == "konstantin"
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "donkon" ## USER on main machine in office
        dt = CSV.read("data/FlowDataGermans.csv", DataFrame)
    elseif ENV["USER"] == "dlakelan"
        dt = CSV.read("data/simulations.csv", DataFrame)
        DataFramesMeta.@transform!(dt,:flows = round.(Int32,:predict),:frompop_ger = :frompop, :topop_ger = :topop)
    end
end    

make_random_inits = function(Nages, ncoefs)
    inits = [
        rand(Normal(0.0, 1.0), Nages); #a
        rand(Gamma(3.0, 1.0 / 2.0), Nages); #b
        rand(Uniform(1.5, 2.5), Nages); #c
        rand(Uniform(0.0, 0.03), Nages); #d0
        [1.5, 8.0]; # neterr and logisticconst
        fill(0.0, Nages); #kd
        fill(0.0, Nages * ncoefs) # desirecoefs
        ]
  return inits
end