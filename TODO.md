# TODO

- I rm Manifest.toml from the repo, because I read people
   instantiating our julia env need everything in very specific
   versions then, is this correct? One issue now is when switching
   branches that have different packages installed, then Manifest and
   Project won't agree
- Share folder with Daniel with all writing output
   - overview
   - paper
   - analysis
- Can we use RBF's coefs to speed up fitting? Check if they are
   correlated at all between years
- If the very dark spots in early years are due to data issues, why is
  there no very bright spot? All the people moving from Göttingen need
  to go somewhere. Did I rm some rows from all the data?
- Not important cause it works, but it triggers a warning. So I added
  @suppress. rm age and year from chain and save them as metadata
  of chain, if possible, would help. If not, just return (; chn, mdl, prd, age,
  year)
- Define new struct OptimResult and declare this type for analyze() and so forth
- in helpeR there are scripts not actually used, like utils.R, plots.R
  ..., rm and clean repo before release
- the cleaning and preprocessing of density data, where did I do this?
- Same with shapes, where did I download and cleaned them?
