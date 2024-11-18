function usdiagplots(alldata,paramvec,parnames)

    bar(1:length(paramvec),paramvec; title="Parameter Values",xlab="index") |> display
    #mapest=deserialize("fitted_models/USmodel_map_2024-08-29T17:24:15.073.dat")
    (densmin,densmax) = (alldata.model.args.densmin, alldata.model.args.densmax)
    (xmin,xmax) = (alldata.model.args.xmin, alldata.model.args.xmax)
    (ymin,ymax) = (alldata.model.args.ymin, alldata.model.args.ymax)

    kdindx = (1:36) .+ 6
    desindx = (1:36) .+ (6+36)
    kdfun = Fun(ApproxFun.Chebyshev(densmin .. densmax) * ApproxFun.Chebyshev(densmin .. densmax),paramvec[kdindx] ./ 10)
    desirfun = Fun(ApproxFun.Chebyshev(xmin .. xmax) * ApproxFun.Chebyshev(ymin .. ymax ),paramvec[desindx] ./ 10)

    plotlyjs()# set the background

    p1=heatmap(kdfun;  title="dens Fun (dens,dens)",bins=(40,40),size=(600,600),xlim=(-7,7),ylim=(-7,7)) 

    p2=StatsPlots.histogram2d(StatsBase.sample(alldata.model.args.logreldens,1000),StatsBase.sample(alldata.model.args.logreldens,1000),normalize=true,
        size=(600,600),xlim=(-7,7),ylim=(-7,7),alpha=.25, title="Distribution of logreldens pairs")
    StatsPlots.plot(p1,p2; layout=(2,1)) |> display

    PlotlyJS.plot(PlotlyJS.histogram2dcontour(x=StatsBase.sample(alldata.model.args.logreldens,1000), y=StatsBase.sample(alldata.model.args.logreldens,1000)),
        PlotlyJS.Layout(alpha=0.1, title="Distribution of Randomly chosen Density Pairs",xlab="logreldensity",ylab="logreldensity")) |> display

    heatmap(desirfun, title="Desirability Fun (long,lat)",c=:rainbow) |> display

    preds = generated_quantities(alldata.model,paramvec,parnames)
    alldata.flows.preds = preds

    netactual = usnetmig(levelcode.(alldata.flows.fromcounty),levelcode.(alldata.flows.tocounty),alldata.flows.COUNT)
    netpred = usnetmig(levelcode.(alldata.flows.fromcounty),levelcode.(alldata.flows.tocounty),alldata.flows.preds)

    density(netactual-netpred; title="Residual Net error (actual - pred)",xlab="Count of people") |> display
    density((netactual-netpred) ./ alldata.geog.POPESTIMATE2016; xlim=(-.3,.3),legend=false, title="Residual Net error (actual - pred)",xlab="Fraction of source county") |> display
    density(netactual ./ alldata.geog.POPESTIMATE2016; xlim=(-.3,.3),legend=false,title="Net migration as fraction of population") |> display
    density(netactual ./ 1000.0; title="Actual Net Migration (thousands)",xlim=(-20,20)) |> display
    density(netactual ./ 1000.0; title="Actual Net Migration (thousands, full range") |> display

    scatter(netpred ./ 1000, netactual ./ 1000; ylab="actual net migration (thousands)",xlab="predicted (thousands)",xlim=(-20,20),ylim=(-20,20),
        title= "Net Migration prediction")
    Plots.abline!(1,0; legend=false) |> display

    scatter(netpred ./ 1000, netactual ./ 1000; ylab="actual net migration (thousands)",xlab="predicted (thousands)",
        title = "Net Migration prediction")
    Plots.abline!(1,0) |> display

    samps = StatsBase.sample(eachindex(alldata.flows.dist),5000; replace=false)

    flowsamp = alldata.flows[samps,:]

    density(log.(flowsamp.preds); label = "log(preds)", title="Density of log predictions")    
    density!(log.(flowsamp.COUNT .+ .01); label="actuals + .01") |> display
    StatsPlots.scatter(flowsamp.dist,log.((flowsamp.COUNT .+ .01) ./ flowsamp.preds); alpha=0.1, title="log(flow/pred) vs dist (in km)") |> display

    scatter(log.(flowsamp.preds/median(flowsamp.preds)),log.(flowsamp.COUNT/median(flowsamp.preds)); title="Comparing prediction to actual",
        xlab="log(pred/median(pred))", ylab="log(COUNT/median(pred))",xlim=(-5,5),ylim=(-5,5), alpha=0.1) |> display


    StatsPlots.scatter(flowsamp.dist,[log.((r.COUNT .+ 1) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((flow+1)/from_pop) vs distance", alpha=0.1) |> display

    StatsPlots.scatter(flowsamp.dist,[log.((r.preds) ./alldata.geog.POPESTIMATE2016[levelcode(r.fromcounty)]) for r in eachrow(flowsamp)];
        title = "log((pred+1)/from_pop) vs distance", alpha=0.1) |> display

    kfrom = paramvec[6] / 10.0
    density(kfrom .* log.(alldata.geog.POPESTIMATE2016[levelcode.(flowsamp.tocounty)] ./ median(alldata.geog.POPESTIMATE2016)); title="DIstribution of log population US * kfrom") |> display

    tots = totalflow(alldata.flows.fromcounty,alldata.flows.tocounty,alldata.flows.COUNT,alldata.geog.countyid)
    totpreds = totalflow(alldata.flows.fromcounty,alldata.flows.tocounty,alldata.flows.preds,alldata.geog.countyid)

    scatter(log.(totpreds ./ alldata.geog.POPESTIMATE2016),log.(tots ./ alldata.geog.POPESTIMATE2016),
        xlab="log(Pred total flux / Pop)", ylab = "log(Actual Total Flux / Pop)", title = "Total Flux comparison")
    Plots.abline!(1,0; label = "y = x", legend = false) |> display

end
