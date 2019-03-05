# portalcasting (development version)

# [portalcasting 0.3.1](https://github.com/weecology/portalcasting)
*active development*


# [portalcasting 0.3.0](https://github.com/weecology/portalcasting/releases/tag/v0.3.0)
*2019-03-04*

### Completed migration of plotting code
* `plot_cast` is now `plot_cast_ts` and is now fully vetted and tested
* `plotcastts_ylab` and `plotcastts_xaxis` provide tidied functions for
producing the y label and x axis (respectively) for `plot_cast_ts`.
* `plot_cast_point` is now added to replace `plot_species_forecast`.
* `plotcastpoint_yaxis` provides tidied functionality for the y axis of 
`plot_cast_point`.
* `select_most_ab_spp` allows for a simple selection of the most abundant
species from a -cast.
* `plot_err_lead_spp_mods` and `plot_cov_RMSE_mod_spp` now added to 
replace the raw code in the evaluation page.

### Processing of forecasts
* `read_casts` (old) is now `read_cast` and specifically works for only one -cast.
* `read_casts` (new) reads in multiple -casts.
* `select_cast` is now `select_casts` and allows a more flexible selection
by default.
* `make_ensemble` now returns a set of predictions with non-`NA` bounds when 
only one model is included (it returns that model as the ensemble).
* `most_recent_cast` returns the date of the most recent -cast. Can be dependent
on the presence of a census.
* `verify_cast` and `cast_is_valid` replace `forecast_is_valid` from the 
repo codebase. `verify_cast` is a logical wrapper on `cast_is_valid` that 
facilitates a pipeline integration. `cast_is_valid` does the major set of
checks of the cast data frame.  
* `append_observed_to_cast` is provided to add the observed data to the forecasts
and add columns for the raw error, in-forecast-window, and lead time as well.
* `measure_cast_error` allows for summarization of errors at the -cast level.

### Processing of data
* `most_recent_census` returns the date of the most recent census.

### Minor changes
* Argument order in `models` is reversed (`add` then `set`) and defaults in general
are now `NULL` and `NULL`, but `set = "prefab"` within the options functions, to
make it easy to run a novel model set.
* Argument order in `subdirs` is reversed (`subs` then `type`) and defaults in 
general are now `NULL` and `NULL`, but `type = "portalcasting"` within options
functions and `dirtree` to make it easier to manage a single subdirectory.
* `fdate` argument has been replaced throughout with `cast_date` for generality.

### Utilities
* `na_conformer` provides tidy functionality for converting non-character `NA` 
entries (can get read in from the data due to the `"NA"` species) to `"NA"`. 
Works for both vectors and data frames.

# [portalcasting 0.2.2](https://github.com/weecology/portalcasting/pull/82)
*2019-02-12* 

### Beginning to migrate plotting code
* `plot_cast` is developed but not yet fully vetted and tested, nor integrated
in the main repository. It will replace `forecast_viz` as a main plotting
function.

### Added `"moons"` to `read_data`
* `read_data`'s options have been expanded to include `"moons"`.
* Not fully implemented everywhere, but now available.

### Bug fix in `interpolate_data`
* `interpolate_data` was using `rodent_spp` in a way that assumed the `"NA"` 
species was coded as `"NA."`, which it wasn't. 
* Expansion of `rodent_spp` to include an `nadot` logical argument, with default
value of `FALSE`.

# [portalcasting 0.2.1](https://github.com/weecology/portalcasting/pull/81)
*2019-02-12* 

### Bug fix in `read_data`
* `read_data` was reading the All rodents file for Controls as well, which caused
the forecasts for the Controls to be duplicated of the All forecasts.
* Simple correction here.

# [portalcasting 0.2.0](https://github.com/weecology/portalcasting/pull/79)
*2019-02-04* 

### Code testing
* All of the code is now tested via [`testthat`](https://github.com/weecology/portalcasting/tree/master/tests).
* Test coverage is tracked via [Codecov](https://codecov.io/gh/weecology/portalcasting).
* The only functionality not covered in testing on Codecov is associated with
`download_predictions()`, which intermittently hangs on Travis. Testing is
available, but requires manual toggling of the `test_location` value to
`"local"` in the relevant test scripts (02-directory and 12-prepare_predictions).

### Enforcement of inputs
* Most of the functions previously did not have any checks on input argument 
classes, sizes, etc.
* Now all functions specifically check each argument's value for validity
and throw specific errors.

### Documentation
* All of the functions have fleshed out documentation that specify argument
requirements, link to each other and externally, and include more information.

### Data classes
* To smooth checking of different data structures, we now define data
objects with classes in addition to their existing (`data.frame` or
`list`) definitions.
* `rodents`, `rodents_list`, `covariates`, `metadata`, `moons`
* These classes do not presently have any specified methods or functions.

### Options list classes 
* To smooth checking of different list structures, we now define the options
list objects with classes in addition to their existing `list` definitions.
* `all_options`, `dir_options`, `PortalData_options`, `data_options`, 
`moons_options`, `rodents_options`, `covariates_options`, `metadata_options`,
`predictions_options`, `models_options`, and `cast_options`
* Each of these classes is created by a function of that name.
* These classes do not presently have any specified methods or functions 
that operate on them.

### Added functions
* `classy()` allows for easy application of classes to objects in a `%>%` pipeline
* `read_data()` provides simple interface for loading and applying classes to
model-ready data files.
* `remove_incompletes()` removes any incomplete entries in a table, as defined
by an `NA` in a specific column.
* `check_options_args()` provides a tidy way to check the input arguments (wrt
class, length, numeric limitations, etc.) to the options functions.

### Vignettes
* Three vignettes were added:
  * [current models vignette](https://weecology.github.io/portalcasting/articles/models.html) was brought from the [forecasting website](https://portal.naturecast.org/models.html).
  * [codebase vignette](https://weecology.github.io/portalcasting/articles/codebase.html) was created from the earlier `README.md` file.
  * [adding a model vignette](https://weecology.github.io/portalcasting/articles/adding_a_model.html) was constructed based on two pages from the 
    [Portal Predictions repository](https://github.com/weecology/portalPredictions) wiki
    ([1](https://github.com/weecology/portalPredictions/wiki/Adding-a-new-model) and 
    [2](https://github.com/weecology/portalPredictions/wiki/Forecast-file-format)) and with substantial additional text added.

### Retention of all forecasts of covariates
* Previous versions retained only one covariate forecast per newmoon.
* We now enable retension of multiple covariate forecasts per newmoon and tag 
the forecast date with a time stamp as well.

### Website
* Added a website driven by [pkgdown](https://pkgdown.r-lib.org/).

### Changelog
* Developed this changelog as part of the package.

### Support documents
* Added [code of conduct](https://github.com/weecology/portalcasting/blob/master/CODE_OF_CONDUCT.md) 
and [contribution guidelines](https://github.com/weecology/portalcasting/blob/master/CONTRIBUTING.md)
to the repository.

# [portalcasting 0.1.1](https://github.com/weecology/portalcasting/commit/05ed76f76a82f32a5a3120eb7c9ef0dc95bd8ae4)
*2019-01-13*

### Addressing Portal Data download
* Setting default back to the Zenodo link via updated portalr function.
* Updated `fill_PortalData()` and new `PortalData_options()` allow 
for control over download source.

# [portalcasting 0.1.0](https://github.com/weecology/portalcasting/tree/77fb0c8a32de5ff39715e652ce8e5b813ad02ff3) 
*true code edits 2018-12-14, version number updated 2019-01-02*

### Migration from Portal Predictions repository
* Code was brought over from the forecasting repository to be housed in its 
own pacakge.
* Multiple updates to the codebase were included, but intentionally just
"under the hood", meaning little or no change to the output and simplification
of the input.
* A major motivation here was also to facilitate model development, which 
requires being able to set up a local version of the repository to play with
in what we might consider a ["sandbox"](https://en.wikipedia.org/wiki/Sandbox_(software_development)).
This will allow someone to develop and test new forecasting models in a space
that isn't the forecasting repo itself (or a clone or a fork), but a truly novel
location. At this point, the sandbox setup isn't fully robust from within this 
package, but rather requires some additional steps (to be documented).

### Development of code pipeline
* The previous implementation of the R codebase driving the forecasting 
(housed within the [portalPredictions](https://github.com/weecology/portalPredictions/tree/ac032f3938a6695a8e5d27ee380032195b2af396)
repo) was a mix of functions and loose code and hard-coded to the repo.
* The package implementation generalizes the functionality and organizes
the code into a set of hierarchical functions that drive creation and use
of the code within the repo or elsewhere. 
* See the [codebase vignette](https://weecology.github.io/portalcasting/articles/codebase.html) for further details.

### Explicit directory tree
* To facilitate portability of the package (a necessity for smooth sandboxing
and the development of new models), we now include explicit, controllable
definition of the forecasting directory tree.
* See the [codebase vignette](https://weecology.github.io/portalcasting/articles/codebase.html) for further details.

### Introduction of options lists
* To facilitate simple control via defaults argument inputs and flexibility to 
changes to inputs throughout the code hierarchy, we include a set of functions
that default options for all aspects of the codebase.
* See the [codebase vignette](https://weecology.github.io/portalcasting/articles/codebase.html) for further details.

# [portalcasting 0.0.1](https://github.com/weecology/portalPredictions/tree/ac032f3938a6695a8e5d27ee380032195b2af396) 
*2018-11-06*

* This is the last iteration of the code that now exists in portalcasting in
its previous home within the portalPredictions repo. 
* It was not referred to by the name portalcasting at the time.
