# Project Notes and Issues

## Fitting

Fitting model3 is crazy hard. We've had some issues related to "bugs"
but also it's a nonlinear model and good inits are annoyingly hard to
find. We haven't had as much success dealing with Pluto or interactive
graphic methods, this needs to get fixed so we can find good initial
conditions.

We actually have had some great fits in the past, but it was through a
process of extremely careful constraining chebyshevs and kd values and
narrow windows on certain parameters etc and then widening those
windows a little at a time. This method corresponds to something I
tried to program into the testmod3 function but has been untested.

The existence of those good fits makes me think we can do well with
our model if we were to find the proper region.

A big problem with finding the proper region is, even if we do create
some good initial values, the 800,000 data points which is ~100k data
points per age group leaves us with posterior distributions on things
like a which are likely to be

## reduced data

It's definitely going to be easier to find the optimum region if there
is far less data constraining the problem and so the region is
bigger. This is like taking a bullseye on a target that's a couple mm
wide and blowing it up to at least a few cm. To this end, it would be
really a good idea to fit the model on a selection of perhaps 400
flows total. I suggest splitting into 4 regions of distance, and
subsetting by distance..

100 flows from distances less than 100 km
100 flows from 100-200km
100 flows from 200-400km
100 flows from 400+km

maybe do this for each age group, so total flows is 400*6 = 2400 rows
in the data set. If we can't fit this data, we are well and truly
screwed.
