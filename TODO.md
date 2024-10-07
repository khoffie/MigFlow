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
