Some to-do thoughts:

1) hexbin(fromdens,todens) for zero flows showed us that a lot of them come from and go to low density regions 
We should try to make the probability of a flow drop a bunch when density in the two regions is low. 

2) drop the mixture model and actually model the individual flows using density

3) Create a Julia diagnostic plotting routine, that shows log(flow/pred) vs distance, hexbin(log(fromdens),log(todens)), maps of germany, other stuff

4) Further investigation of why optimization doesn't work well on simulated data.

5) A new simulator for data? Use Turing prior predictive?

6) A Radial Basis Function expansion instead of Chebyshev?



## Bugs and what to improve
1) main() runs even if US data not downloaded
2) needed to create manuscript_input manually
3) write downloadUSdata() that checks if individual data exists and if not downloads
5) various loadUS functions should check if data exists

## Wokflow US and German data
download US data:
	- getUSflows
	- getUSgeog
	- getUScountypop

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
