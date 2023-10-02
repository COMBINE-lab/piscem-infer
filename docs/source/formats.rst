Formats
=======

Quantification output
---------------------

The quantification results are provided in the format of a tab-separated values (TSV) file. The file will 
be named based on whatever output stem is provided to the ``piscem-infer`` command via the ``-o`` flag. 
Specifically, with and output name of ``res``, the quantification results will be written to a file named 
``res.quant``.

The file contains a single header row, and consists of 4 columns with the column names and descriptions as 
specified below:

+----------------------------------------+-----------------------------------------------+--------------------------------------------------------------+---------------------------------------------------------+
| target_name                            | eeln                                          |  tpm                                                         | ecount                                                  | 
+========================================+===============================================+==============================================================+=========================================================+
| The identifer of the quantified target | The expected *effective length* of the target | The abundance of the target in transcripts per million (TPM) | The expected number of fragments assigned to the target |
+----------------------------------------+-----------------------------------------------+--------------------------------------------------------------+---------------------------------------------------------+


Fragment length distribution
----------------------------

Though not directly intended for use by end-users, the fragment length distribution obtained from the 
``RAD`` file used by ``piscem-infer``, and used to compute the expected effective lengths (eelen) of each 
target, is output in a `Apache Parquet <https://parquet.apache.org/>`_ format file named ``res.fld.pq``, where 
``res`` is the output stem provided to the ``-o`` option of ``piscem-infer``.

The Parquet file can be loaded using any of the commonly-available libraries for reading this common format in 
languages like Python and R.


Inferential replicates
----------------------

If ``piscem-infer`` was run with a ``--num-bootstraps`` value greater than 0, then a file called ``res.infreps.pq``
will be written (again, where ``res`` is the output stem provided to the ``-o`` option of ``piscem-infer``). The 
``infreps`` file is also a Parquet format file, and it stores in each *column* the value of a single inferential 
replicate (i.e. an abundance for each target).  The number of columns is equal to the number of requested bootstraps,
and the number of rows is equal to the number of targets represented in the header of the ``RAD`` file.  These 
inferential replicates can be used to assess the inferential uncertainty in the provided quantification estimates.
In other words, variance across a row represents uncertainty in the estimated abundance of the corresponding target.
While the ability to load up this data frame easily in Python or R makes this information readily available to 
end-users, its primary purpose is to be used in uncertainy-aware methods for differential analysis (such as 
`swish <https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html>`_, to which we hope 
to add ``piscem-infer`` support soon).


.. autosummary::
   :toctree: generated

    piscem-infer
