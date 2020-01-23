# Changelog

## Unrealeased

### Added
- Added support for BED files without comma terminator in fields 11 and 12
- Added --transcript_feature_name argument for gtf2bed

## [0.2.2] - 2019/02/27

### Added
- Better documentation and examples of API methods

### Fixed
- Fixed shadowing of built-in `all` function by translateChr()

## [0.2.1] - 2019/02/25

### Added
- Added support for strandedness in tx2genome()

### Fixed
- Fixed bug in tx2genome() strand handling

## [0.2.0] - 2019/01/19

### Fixed
- Improvements to the CLI messages
- Major documentation rewrite

## [0.1.2] - 2018/11/09

### Added
- Added --filterKey option to gtf2bed
### Fixed
- Fixed exception handling in tx2genome
- Corrected off-by-one error in tx2genome
- Added unit tests for tx2genome
- Various documentation fixes

## [0.1.1] - 2018/11/12

### Added
- Added --whichExon option to bed12tobed6

### Fixed
- Fixed a bug in bed12tobed6


