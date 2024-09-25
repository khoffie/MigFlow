## Output Structure I need
I think three tables will be enough. One table for origin-destination
data, like flow predictions. Another table for region data (not
paired) like net prediction and attractiveness. And the last for
parameter estimates.

How should they be structured?

### Origin-Destination Table
Columns: fromdist, todist, year, agegroup, distance, actual, preds

I recommend adding a year column even though we did not use it yet. It
will avoid confusion if we later extent the model to other years and
also a year column in in the German data. I need it to join additional
data like gdp and stuff

### Region Table
Columns: distcode, year, agegroup, actualnet, prednet, actualtotal,
predtotal, cheby_density, cheby_geography

If calculating actualnet, prednet, actualtotal and predtotal is
inconvenient for you, no worries, I can do it easily myself. I am
merely assuming that you generate those anyways.

### Parameter estimates
Columns: agegroup, param_name, estimate

I think this should be it. I really only need flow predictions, the
two cheby predictions, unique identifiers wrt. to age, year and
district / district pair and the parameter estimates per agegroup. The
other stuff I can calculate myself. See what works for you and thank
you!
