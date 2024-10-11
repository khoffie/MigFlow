Some to-do thoughts:

1) hexbin(fromdens,todens) for zero flows showed us that a lot of them come from and go to low density regions 
We should try to make the probability of a flow drop a bunch when density in the two regions is low. 

2) drop the mixture model and actually model the individual flows using density

3) Create a Julia diagnostic plotting routine, that shows log(flow/pred) vs distance, hexbin(log(fromdens),log(todens)), maps of germany, other stuff

4) Further investigation of why optimization doesn't work well on simulated data.

5) A new simulator for data? Use Turing prior predictive?

6) A Radial Basis Function expansion instead of Chebyshev?

# Paper
So, we developed a model that does a lot differently than current
approaches. What are the two most important goals of our writing?

1. Describing clearly what we did and why we did it.
2. Describing clearly how this relates to standard approaches.

Many papers only do 1. and neglect 2. This is a mistake, if you do new
and complicated stuff, you have to explain to readers why they should
care. This responsibility lies with the authors, not the readers. They
have to understand and understanding is easier if someone talks about
things you already know.

Other than that, our paper should be as simple as possible. I suggest
the following:

- We write a methods paper. We focus on the model, show how it works
and that interesting things follow, like two populations: distance
sensitive and insensitive and decoupling of density and geographic
attractiveness.
- We focus on three or so key plot types: distance dependence,
density hex plot, attractiveness maps. 
- We do not relate our findings closely to internal migration studies,
it would be too much. Same with careful interpretations of the model's
implications (like this change in attractiveness is related to rising
rents, whereas here a factory closed). I guess this is clear if we
write a methods paper.

## Outline of Paper
1. Introduction: What do we write here? We could focus on:
   - Debates around gravity model: Do not explain anything v are
     empirically very successful. Always good to start with an
     apparent paradox!
   - Concrete improvement of gravity models: Distance sensitive v
     distance insensitive
   - Scientific v Statistical model: Nobody has any idea how to think
     scientifically about functional relationship. We do! Maybe some
     abstract thoughts about science in general
2. Main Part / Methods
   - Comparing our d0 model with gravity model. Explain new distance
     dependence, plot log(actual / pred) v distance
   - Adding forth gravity variable density. Explaining, density hex plots
   - Geographical attractiveness
3. Summary / Discussion
   - ?

## Who writes what?
So, how do we write the paper? I like to write and I have to write
about our model anyways for my thesis. However, the model is based
mainly on your thoughts. You are likely much better at explaining the
model. And you are a native speaker. What do you think? I could write
about German context and since I know the literature better about the
context.

## Open Questions

## Further Thoughts
While this seems straightforward, not everything is clear.

1. I presented the gravity model as standard approach and I think it
   still is, but clear trends are visibly going in a different (and
   entirely wrong, if you ask me) direction. Econometricians like this
   bivariate flow data because "it allows a rich structure for fixed
   effects". Meaning, they often add origin and destination fixed
   effects and sometimes even for origin-destination-pairs. Then, they
   claim identification of causal effects is better. Does not seem to
   further understanding of the processes involved.
2. How much attention do we pay to comparing US and Germany? In
   general I like the idea. It is interesting to see how inferences
   differ, but doing it well would be a new paper I feel. There are
   some issues: In Germany we have agegroups, in the US we only fit
   positive flows. 
3. Gravity models are used in internal migration modeling,
   international migration modeling and modeling of international
   trade flows. Our critique of the gravity model is partly general:
   Coefs on pop, singularity in distance = zero. Density and
   geographic ideas are more specific to internal migration, so I
   think we choose this as topic
4. We always speak about internal migration. The literature is not
   clear on what this means. It is also common, and useful, to
   differentiate between residential mobility (short distance moves)
   and internal migration (long distance). Our model is really about
   both, we only have to stick to one term and mention it refers to
   moves across all distances.

# Workflow
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
