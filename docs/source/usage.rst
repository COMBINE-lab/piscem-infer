Usage
=====

.. _installation:

Installation of ``piscem-infer``
--------------------------------

To install piscem-infer, use `cargo <https://github.com/rust-lang/cargo>`_:

.. code-block:: console

   $ cargo install piscem-infer


or install it from source


.. code-block:: console

   $ git clone https://github.com/COMBINE-lab/piscem-infer
   $ cd piscem-infer
   $ cargo build --release



Command line parameters for ``piscem-infer``
--------------------------------------------

Currently, ``piscem-infer`` has only a ``help`` and ``quant`` sub-command. The ``help`` sub-command is 
self-descriptive, and the usage for the ``quant`` sub-command is provided below.

.. code-block:: console

   quantify from the rad file

   Usage: piscem-infer quant [OPTIONS] --input <INPUT> --lib-type <LIB_TYPE> --output <OUTPUT>

   Options:
   -i, --input <INPUT>                              input stem (i.e. without the .rad suffix)
   -l, --lib-type <LIB_TYPE>                        the expected library type
   -o, --output <OUTPUT>                            output file prefix (multiple output files may be created, the main will have a `.quant` suffix)
   -m, --max-iter <MAX_ITER>                        max iterations to run the EM [default: 1500]
         --convergence-thresh <CONVERGENCE_THRESH>    convergence threshold for EM [default: 0.001]
         --presence-thresh <PRESENCE_THRESH>          presence threshold for EM [default: 0.00000001]
         --param-est-frags <PARAM_EST_FRAGS>          number of (unique) mappings to use to perform initial coarse-grained estimation of the fragment length distribution. These fragments will have to be read from
                                                      the file and interrogated twice [default: 500000]
         --factorized-eqc-bins <FACTORIZED_EQC_BINS>  number of probability bins to use in RangeFactorized equivalence classes. If this value is set to 1, then basic equivalence classes are used [default: 64]
         --fld-mean <FLD_MEAN>                        mean of fragment length distribution mean (required, and used, only in the case of unpaired fragments)
         --fld-sd <FLD_SD>                            mean of fragment length distribution standard deviation (required, and used, only in the case of unpaired fragments)
         --num-bootstraps <NUM_BOOTSTRAPS>            number of bootstrap replicates to perform [default: 0]
         --num-threads <NUM_THREADS>                  number of threads to use (used during the EM and for bootstrapping) [default: 16]
   -h, --help                                       Print help
   -V, --version                                    Print version


Most of the parameters are self-explanatory.  The ``--input`` option should point to the stem of the input ``RAD`` file, the ``output`` should point to the output stem.  This can contain directories (which
will be created if they do not yet exist, and several files with this prefix stem, but different suffixes, will be created). 

The ``--lib-type`` option describes the library type (how the reads are expected to have arisen from the underlying molecules), and is specified 
using `salmon's library type specification <https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype>`_. The ``--param-est-frags`` determines how many fragments are used from the input
RAD file to estimate the empirical fragment length distribution.  Not too many fragemnets are necessary to estimate the FLD reliably, but these fragments will be interrogated twice, so this parameter can have 
a (generally small) effect on the runtime.  The ``--factorized-eqc-bins`` option determines how many conditional probability bins are used in the range factorized equivalence classes used for inference.
To learn more about range factorized equivalence classes and the choice of this parameter, please see the paper 
`Improved data-driven likelihood factorizations for transcript abundance estimation <https://doi.org/10.1093/bioinformatics/btx262>` which describes the purpose of this improved notion of equivalence 
classes and the effect of choosing more or fewer bins.  If the number of bins is set to 1, then there is no purpose in using range factorized equivalence classes, and instead basic equivalence classes
are used.  If the sample being processed is not paried-end, then an empirical FLD cannot be reliably estimated from the data, and the ``--fld-mean`` and ``--fld-sd`` should be provided by the user, which 
will be used as the mean and standard deviation of a truncated normal distribution of fragment lengths.  Setting ``--num-bootstraps`` to some value ``k > 0`` will cause ``k`` inferential replicates
to be generated.  Finally ``--num-threads`` specifies the number of threads that will be used during inference (to run the main estimation EM algorithm as well as to perform bootstrapping).
Currently, the parsing of the RAD file is single-threaded, and so these threads are only used during inference an bootstrapping.



An example run
--------------

Here, we demonstrate how to process some data, that can then be used for a simple differential 
expression analysis.  We'll keep all of our analysis in a single directory that we'll create now.

.. code-block:: console

   $ mkdir piscem_tutorial
   $ cd piscem_tutorial


Obtaining the reference
~~~~~~~~~~~~~~~~~~~~~~~

First, we'll grab the reference transcriptome we will use for quantification - in this case Gencode's v44 annotation 
of the human transcriptome.

.. code-block:: console

   $ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

Building the index
~~~~~~~~~~~~~~~~~~

Next, we'll build the index. You'll only have to do this once (or whenever you want to update the annotation you're using). To 
build the index and map the reads, we'll need ``piscem``. You can either build it from source according to the instructions 
on the `GitHub page <https://github.com/COMBINE-lab/piscem>`_, or you can install it from ``biconda`` using ``conda install piscem``. 
Once you have it installed, you can build the index with:

.. code-block:: console

    $ piscem build -s gencode.v44.transcripts.fa.gz -k 31 -m 19 -t 16 -o gencode_v44_idx

Obtaining the reads
~~~~~~~~~~~~~~~~~~~

To obtain some sample read data, we'll use the excellent |fastqdl|_ tool that you can install 
via either ``pip`` or bioconda (through ``conda`` or ``mamba``).

.. code-block:: console
    
   $ accessions=(SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521)
   $ for sra in ${accessions[@]}; do fastq-dl -a $sra; done

This will retrieve 8 accessions (16 files in total since each sample is a paired-end sequencing run).

Mapping the reads
~~~~~~~~~~~~~~~~~

Next, we'll use ``piscem`` again to map the reads.  The following command will do it for us (you can check out ``piscem map-bulk -h`` for 
a descripton of all the options):

.. code-block:: console

   $ mkdir -mappings
   $ for acc in ${accessions[@]}; do 
      piscem map-bulk -t 16 -i gencode_v44_idx -1 ${acc}_1.fastq.gz -2 ${acc}_1.fastq.gz -o mappings/${acc}


Quantification with ``piscem-infer``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we've mapped the reads to produce a bulk-oriented ``RAD`` file, we're ready to quantify with ``piscem-infer``!
Here, in addition to performing the basic quantification, we will be creating inferential replicates (i.e. bootstrap
samples) for each sample we quantify. This is designated by the ``--num-bootstraps`` parameter. To perform the 
bootsrapping in parallel, we'll make use of multiple threads (``--num-threads 16``).

.. code-block:: console
  
   $ for acc in ${accessions[@]}; do
      piscem-infer quant --num-bootstraps 16 --num-threads 16 -i mappings/${acc} -l IU -o quant/${acc}

Note that we pass to the ``-o`` flag a file *stem* prefixed with a path (in this case ``quant``). This is because ``piscem-infer``
will produce several output files.  All of them will share the same *stem*.  If we pass a stem that is prefixed with some path 
(e.g. a directory) then this directory will be created if it doesn't exist. We also let ``piscem-infer`` know the library type 
(i.e. how we expect the reads to map), where ``piscem-infer`` uses `salmon's library type specification <https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype>`_.
Here we expect the library to be unstranded and the paired-end reads to map "inward" (i.e. facing each other).

If we look at the files generated with the stem corresponding to, say, the second sample (``SRR1039509``), we 
see the following:

.. code-block:: console

    $ ls -la quant/${accessions[1]}*
      .rw-rw-r--@ 3.1k rob  5 Oct 15:12 quant/SRR1039509.fld.pq
      .rw-rw-r--@  12M rob  5 Oct 15:15 quant/SRR1039509.infreps.pq
      .rw-rw-r--@ 1.1k rob  5 Oct 15:15 quant/SRR1039509.meta_info.json
      .rw-rw-r--@  36M rob  5 Oct 15:12 quant/SRR1039509.quant 

The file ``SRR1039509.quant`` contains the quantification estimates, and is of a very similar format to e.g. a ``salmon`` ("quant.sf") format file.  The file format for the quantification result, as well as that of other outputs, is described in the :ref:`format section of this documentation<Quantification output>`. The file ``SRR1039509.meta_info.json`` contains 
information about the quantification run.  The files ``SRR1039509.fld.pq`` and ``SRR1039509.infreps.pq`` are `Apache Parquet <https://parquet.apache.org/>`_ format files and contain, respectively, information about the inferred fragment length distribution of the sample and the inferential replicates that we requested to be computed.


Subsequent differential analysis using ``tximport`` and ``DESeq2``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we'll show how to perform differential analysis (at the gene level) with the quantification 
estimates we just computed using ``tximport`` and ``DESeq2``.  First, we'll need 
*just a little bit more information*. We'll need a file containing the run information about these 
samples (which includes, e.g. the metadata about how they were treated), and a file containing the 
transcript-to-gene mapping. To make this tutorial easier to follow, these can be obtained directly 
using the following commands (we'll download them into our current working directory, where we will 
also perform our differential analysis).

.. code-block:: console

  $wget -O SraRunTable.txt -r --no-check-certificate 'https://drive.google.com/uc?export=download&id=1Qt93SG0rAI-GJ9LCmyl1gqM-9T5JlJBx'
  $wget -O t2g.csv -r --no-check-certificate 'https://drive.google.com/uc?export=download&id=1fUpx-0HHI8msRaZm2UUKf-d5lDD0gYXZ'

Now, we're ready to perform our DE analysis. That part of the tutorial can be found 
in this `Quarto document <_static/simple_de_example.html>`_.


.. |fastqdl| replace:: ``fastq-dl``
.. _fastqdl: https://github.com/rpetit3/fastq-dl

