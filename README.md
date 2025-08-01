# piscem-infer

### What is `piscem-infer`?

`piscem-infer` is a tool to consume bulk-RAD files (produced by [`piscem`](https://github.com/COMBINE-lab/piscem) or [`piscem-cpp`](https://github.com/COMBINE-lab/piscem-cpp)) 
and to produce from them abundance estimates of the targets against which the reads were mapped.  For example, a basic RNA-seq pipeline could consist of mapping the reads 
against the transcriptome using `piscem`, and then quantifying transcript abundance using `piscem-infer`.  Likewise, one could use the pair of tools on metagenomic index and 
metagenomic sequencing reads to perform taxonomic abundance estimation.  The main goal of `piscem-infer` is to separate the statistical inference algorithms and code from the 
code that performs indexing and mapping.  This allows faster and more agile development of new improvements to the inference method, as well as eases the maintenance burden.

The `piscem-infer` program is written in `rust`, which makes it easy for end-users to compile, and which also makes it easy for us to deploy without the need for end-users 
to compile it (it's easy to create statically-linked, pre-compiled executables, and to put the package on `bioconda`). At the same time, this gives us access to the safety 
guarantees of `rust`, making the code easier to develop and maintain with confidence while retaining the efficiency of a high-performance, statically-typed, compiled language.

While `piscem-infer` is in *active development*, it is already very usable! We encourage folks who are curious to try it out, to open GitHub issues for any questions or feature
requests, and to provide feedback on how you'd like to see this next "evolution" of (bulk sequencing) abundance estimation tools evolve!


### Check out the documentation

Detailed documentation for `piscem-infer` is available in our dedicated [ReadTheDocs page](https://piscem-infer.readthedocs.io/en/latest/index.html)
