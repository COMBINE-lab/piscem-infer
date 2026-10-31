# Changelog

All notable changes to this project are documented here. The format follows
[conventional commits](https://www.conventionalcommits.org); entries are
generated from commit messages by [git-cliff](https://git-cliff.org).

## [0.8.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.7.0...v0.8.0) (2026-10-31)


### Bug Fixes

* Bump libradicl 0.10 -> 0.21 to read RAD v2 files ([95a69e3](https://github.com/COMBINE-lab/piscem-infer/commit/95a69e39729693945b52577c63ef1deb658e6c2e))
* Default --em-accel to none (match prior releases and salmon) ([f47c115](https://github.com/COMBINE-lab/piscem-infer/commit/f47c115481f84c8d714f2447c5613b494f076159))
* Coalesce orientation-collapsed classes so packed order is deterministic ([e3862c6](https://github.com/COMBINE-lab/piscem-infer/commit/e3862c6d0a1cc8ef969320348a7107049a1913d2))
* Skip target-less classes from strand-filtered fragments ([5a055f1](https://github.com/COMBINE-lab/piscem-infer/commit/5a055f1cd22609bf07580b9c9fd541a188f9145f))
* Use usize CSR offsets in PackedEqMap; size CSR arrays exactly ([792802c](https://github.com/COMBINE-lab/piscem-infer/commit/792802c8a662674d6484e3efca5a24ca1610c73c))
* Correct --fld-mean/--fld-sd help and log message typos ([e761343](https://github.com/COMBINE-lab/piscem-infer/commit/e7613430a6eac0f3ae5416f89a2e5deca1d85bb2))
* Accept -q after the subcommand, add -t, reject --num-threads 0 ([2732bf8](https://github.com/COMBINE-lab/piscem-infer/commit/2732bf8e93060cb0dcb8263890eb9381ed7a1314))
* Resolve -i the way piscem map-bulk names its outputs ([62b0f31](https://github.com/COMBINE-lab/piscem-infer/commit/62b0f3151693a9a49811a9944f3116685c7c023a))
* Validate the bulk RAD layout and file tags before reading ([e5ba41b](https://github.com/COMBINE-lab/piscem-infer/commit/e5ba41b553f4326e4ca9309b6501e25b548f38f6))
* Report truncated and unfinished RAD files instead of aborting ([8b7cff4](https://github.com/COMBINE-lab/piscem-infer/commit/8b7cff4120f3ee12d91996535f49355d3f2785f4))
* Drop a stray quote from the paired-end FLD warning ([df19d02](https://github.com/COMBINE-lab/piscem-infer/commit/df19d0295b4f6b11f06455aa6d6cb2d351c3b9cc))
* Don't index past the conditional means for a 65,536-nt reference ([f0a9b47](https://github.com/COMBINE-lab/piscem-infer/commit/f0a9b47b55507f7bd94c5ccf5a3d9cce8508b95a))
* Keep the parametric cdf flat past its 65,536-nt table ([63f1027](https://github.com/COMBINE-lab/piscem-infer/commit/63f1027e41e3378a70df8c50854ebaa89f093f35))
* Replace unmaintained arrow2 with parquet 60.1; clear audit findings ([c934915](https://github.com/COMBINE-lab/piscem-infer/commit/c9349152c9702794bcb4f0b2d3b8c01cb5293cf5))

### Build & CI

* Build, test, lint and MSRV workflow; declare rust-version 1.91 ([1d36f27](https://github.com/COMBINE-lab/piscem-infer/commit/1d36f271305389a8927b5c17903c74d1a28517f5))
* Bump_and_publish.sh + git-cliff instead of release-please ([cb0575b](https://github.com/COMBINE-lab/piscem-infer/commit/cb0575bcdc0c0dd4b89fac7ab8ee97da20e1046a))
* Depend on the published libradicl 0.21.2 ([0815846](https://github.com/COMBINE-lab/piscem-infer/commit/081584666b09d60cd6587ce9af6f3ab52cd70967))
* Salmon-infer 2.9.1 ([d5872e6](https://github.com/COMBINE-lab/piscem-infer/commit/d5872e667c88c49dba93df3be271a716190979df))

### Documentation

* Don't claim --em-accel schemes reach the same fixpoint faster ([119de19](https://github.com/COMBINE-lab/piscem-infer/commit/119de19c4a2a3847655a8cfa39e588df9437163b))
* Fix tutorial commands, .quant columns and RAD v2 compatibility ([9dcf118](https://github.com/COMBINE-lab/piscem-infer/commit/9dcf11821bb7d18c7f474704733cfc167e12abf0))

### Features

* Port EM/bootstrap/Gibbs to salmon-infer's deterministic sharded EM ([9feceba](https://github.com/COMBINE-lab/piscem-infer/commit/9fecebabde9d0821668bdbf1ab02c7fc448f2073))
* Adopt salmon's EM convergence criteria ([ff898ce](https://github.com/COMBINE-lab/piscem-infer/commit/ff898ceb9b3f7a545eb512f7ceeb29f617d8d567))
* Warn on EM non-convergence and record it in meta_info.json ([f832cec](https://github.com/COMBINE-lab/piscem-infer/commit/f832cec94e254f0f17bb7bc154a4e744fe7da92a))
* --vbem, salmon's variational Bayesian EM, as an option ([11465c2](https://github.com/COMBINE-lab/piscem-infer/commit/11465c2d901dcc57de6d920110285ee3387ab56d))

### Miscellaneous

* Fmt/clippy cleanup of pre-existing issues (untouched modules) ([4db0769](https://github.com/COMBINE-lab/piscem-infer/commit/4db0769bb42c78e1c4f2449177deeb12299e240e))

### Performance

* Key bulk eq classes without fragment orientations ([5b87fe5](https://github.com/COMBINE-lab/piscem-infer/commit/5b87fe588ffd343004402131de559e8471e28ad6))
* Count known eq classes without allocating a label ([479ae28](https://github.com/COMBINE-lab/piscem-infer/commit/479ae28ef1a1e13593f11e7b415bcc9ca21f3d13))
* Parse the RAD file's chunks on --num-threads threads ([cbd781d](https://github.com/COMBINE-lab/piscem-infer/commit/cbd781da7fe6d84c7780853107141a76b077017f))
* Buffer the .quant writer ([e335092](https://github.com/COMBINE-lab/piscem-infer/commit/e335092c76837811ad68dd94ef6a359764abe260))
* Read bulk RAD records in place; build on salmon-infer 2.9.0 ([e9caa58](https://github.com/COMBINE-lab/piscem-infer/commit/e9caa58ca09421cac9e17e0e85289a359387423e))

### Testing

* Make lib_types_parse fail on parse errors; pin replicate reproducibility ([4de0510](https://github.com/COMBINE-lab/piscem-infer/commit/4de0510443e89cf234d17e8e0311422032143c63))

## [0.6.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.5.2...v0.6.0) (2024-02-26)


### Features

* working using new libradicl ([72120a6](https://github.com/COMBINE-lab/piscem-infer/commit/72120a69dea72c823719495d77e9c2fe3fd934c7))

## [0.5.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.4.0...v0.5.0) (2023-10-10)


### Features

* speed up EM ([7602253](https://github.com/COMBINE-lab/piscem-infer/commit/7602253c656664561b0e03db58c4a94e2d883683))

## [0.4.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.3.0...v0.4.0) (2023-10-05)


### Features

* add length to output file ([6a6a6b0](https://github.com/COMBINE-lab/piscem-infer/commit/6a6a6b06dcc452251aee64e2a2ff413a13e68bc3))

## [0.3.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.2.0...v0.3.0) (2023-10-02)


### Features

* propagate reference provanence signatures through quant if they exist ([890c652](https://github.com/COMBINE-lab/piscem-infer/commit/890c6523e587027a3c4b811d5d0046cf3fff63d1))

## [0.2.0](https://github.com/COMBINE-lab/piscem-infer/compare/v0.1.1...v0.2.0) (2023-10-01)


### Features

* add tpm update docs ([648ed87](https://github.com/COMBINE-lab/piscem-infer/commit/648ed870d72044a9a22ae792099302b4051acedf))


### Bug Fixes

* bump version in tead the docs ([d63d8ba](https://github.com/COMBINE-lab/piscem-infer/commit/d63d8bafaae0c6ba6a5a779a448240d594fe2007))

## [0.1.1](https://github.com/COMBINE-lab/piscem-infer/compare/v0.1.0...v0.1.1) (2023-09-30)


### Bug Fixes

* missing file ([0118da7](https://github.com/COMBINE-lab/piscem-infer/commit/0118da7dc79078d313d783ab9e352f398f9b93ce))

## 0.1.0 (2023-09-30)


### Features

* initial release ([376d3e6](https://github.com/COMBINE-lab/piscem-infer/commit/376d3e650b4cd1a090ee9d64e3e82d4abc1f4aac))
