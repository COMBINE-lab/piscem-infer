Usage
=====

.. _installation:

Installation
------------

To install piscem-infer, use `cargo <https://github.com/rust-lang/cargo>`_:

.. code-block:: console

   $ cargo install piscem-infer

or install it from source

.. code-block:: console

   $ git clone https://github.com/COMBINE-lab/piscem-infer
   $ cd piscem-infer
   $ cargo build --release


An example run
--------------

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

To obtain some sample read data, we'll use the excellent ```fastq-dl`` <https://github.com/rpetit3/fastq-dl>`_ tool that you can install 
via either ``pip`` or bioconda (through ``conda`` or ``mamba``).

.. code-block:: console

   $ fastq-dl -a SRR1039508

Mapping the reads
~~~~~~~~~~~~~~~~~

Next, we'll use ``piscem`` again to map the reads.  The following command will do it for us (you can check out ``piscem map-bulk -h`` for 
a descripton of all the options):

.. code-block:: console

   $ piscem map-bulk -t 16 -i gencode_v44_idx -1 SRR1039508_1.fastq.gz -2 SRR1039508_1.fastq.gz -o SRR1039508_mapped


Quantification with ``piscem-infer``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we've mapped the reads to produce a bulk-oriented ``RAD`` file, we're ready to quantify with ``piscem-infer``!
There exist other options we can pass (e.g. we can perform bootstrap sampling to produce inferential replicates, in which 
case we can also request the use of multiple threads, but here we just invoke the most basic quantification process).

.. code-block:: console

   $ piscem-infer quant -i SRR1039508_map -l IU -o quant/SRR1039508

Note that we pass to the ``-o`` flag a file *stem* prefixed with a path (in this case ``quant``). This is because ``piscem-infer``
will produce several output files.  All of them will share the same *stem*.  If we pass a stem that is prefixed with some path 
(e.g. a directory) then this directory will be created if it doesn't exist. We also let ``piscem-infer`` know the library type 
(i.e. how we expect the reads to map), where ``piscem-infer`` uses `salmon's library type specification <https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype>`_.
Here we expect the library to be unstranded and the paired-end reads to map "inward" (i.e. facing each other).

.. code-block:: console

    $ ls -la quant/
    .rw-rw-r--@ 3.1k rob 30 Sep 12:33 SRR1039508.fld.pq
    .rw-rw-r--@  628 rob 30 Sep 12:33 SRR1039508.meta_info.json
    .rw-rw-r--@  33M rob 30 Sep 12:33 SRR1039508.quant

The file ``SRR1039508.quant`` contains the quantification estimates, and is of a very similar format to e.g. a ``salmon`` ("quant.sf") format file.

