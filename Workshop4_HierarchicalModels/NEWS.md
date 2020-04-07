I changed :

* data/ folder at roots ignored by git for data reading
* chunk names in camel, e.g. `loadData` instead of `load_data`, as it creates bug in cross-referencing
* a bit of text (sorry I rpefer narrative documents)
* true table with `kable()`
* table and legend captions
* data fromatting with `dplyr` and `scalel` in `loadData`
* `save_warmup = F,` to save a bit of R space
* sub-level header to help reading
* table of model speed comparison
