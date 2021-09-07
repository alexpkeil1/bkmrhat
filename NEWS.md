# `bkmrhat` v1.1.1

## Major changes
- N/A

## Minor changes
- Change burnin default to NULL in kmbayes_combine (was 0) to fit typical usage pattern and match documentation
- Added kmbayes_combine_lowmem function, which emulates kmbayes_combine and may work in certain settings where memory constraints result in a "low memory" error in kmbayes_combine.

## Bug fixes
- MINOR: fixed "iters" objects in kmbayes_combine and kmbayes_combine_lowmem output ("bkmrplusfit" objects)

# `bkmrhat` v1.0.2

## Major changes
- N/A

## Minor changes
- N/A

## Bug fixes
- MAJOR: seed change in v1.0.1 introduced a bug that would result in chains having identical seeds. This is fixed. 

# `bkmrhat` v1.0.1

## Major changes
- N/A

## Minor changes
- Removed a warning from the future package about possibly invalid seed values (Set seed=TRUE in underlying `future` function calls, per package instructions)

## Bug fixes
- N/A


# `bkmrhat` v1.0.0

## Major changes
- Added parallel diagnostics and re-packaged bkmr package functions for parallel chains
- Added "continue" functions to continue fitting from an existing bkmr fit

## Bug fixes
- N/A

# `bkmrhat` v0.1.16

## Major changes
- First CRAN release

## Bug fixes
- N/A
