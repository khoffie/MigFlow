# TODO

1. Use maximum_a_posteriori instead of maximum_likelihood
2. Add small constants $\epsilon$ and $\delta$ to model
3. I rm Manifest.toml from the repo, because I read people
   instantiating our julia env need everything in very specific
   versions then, is this correct? One issue now is when switching
   branches that have different packages installed, then Manifest and
   Project won't agree
4. Create a nice makie map   
5. Share folder with Daniel with all writing output
   - overview
   - paper
   - analysis
6. Create and save predictions directly when estimating the model,
   instead of a separate step now
7. Can we use RBF's coefs to speed up fitting? Check if they are
   correlated at all between years
