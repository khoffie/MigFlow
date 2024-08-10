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
