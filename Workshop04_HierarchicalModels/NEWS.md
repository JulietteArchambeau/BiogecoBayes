# Sylvain changes

* data/ folder at roots ignored by git for data reading
* chunk names in camel, e.g. `loadData` instead of `load_data`, as it creates bug in cross-referencing
* a bit of text (sorry I rpefer narrative documents)
* true table with `kable()`
* table and legend captions
* data fromatting with `dplyr` and `scalel` in `loadData`
* `save_warmup = F,` to save a bit of R space
* sub-level header to help reading
* table of model speed comparison
* lognormal model with lognormal priors in fixed effects
* lognormal non centered parameterisation


# Content

**A unique model**

1. Fixed effect 
    1. Distributions
        1. Normal
            1. ...
        1. Normal + log(y)
        1. Lognormal on R
        1. Lognormal on R+
    1. Priors on Normal + log(y)
        1. ...
    1. Vectorisation on Normal + log(y) with best priors
        1. ...
1. One-varying intercept
    1. Normal + log(y) Centered
    1. Normal + log(y) Non Centered
    1. Lognormal on R+ Non Centered
1. Two-varying intercept  
1. Varying intercepts & slopes



