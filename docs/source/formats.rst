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

Meta information about the run
------------------------------

``piscem-infer`` will also collect and output metadata about each run.  This information is written to a file 
``res.meta_info.json``, where ``res`` is the output stem provided to the ``-o`` option of ``piscem-infer``.
The information in this file includes statistics about what was encoded in the input ``RAD`` file, how 
``piscem-infer`` itself was invoked, as well as relevant information about the provenance of the reference 
sequence against which the reads were mapped (if the ``piscem`` index, itself, was built to contain this 
information).  Below, is an illustrative example of the meta information for a single run (the one 
outlined in the :ref:`usage example<An example run>`).  In general, between versions of ``piscem-infer``, new 
fields may be added to the meta-information that is written out, but we intend to be cautious about removing 
or renaming existing fields, since downstream analyses may come to depend on them.  That being said, if there is
information in this file that you are using downstream, or there is other information not in this file that you would 
like to have provided, please let us know.

.. code-block:: json
  {
    "mapped_frag_stats": {
      "filtered_ori_count": [
        0,
        0,
        0,
        0,
        0,
        0,
        0
      ],
      "mapped_ori_count": [
        0,
        255966,
        282047,
        10614009,
        10635715,
        0,
        0
      ],
      "num_mapped_reads": 21325482,
      "tot_mappings": 125259802
    },
    "num_bootstraps": 0,
    "num_targets": 252797,
    "quant_opts": {
      "convergence_thresh": 0.001,
      "fld_mean": null,
      "fld_sd": null,
      "input": "/mnt/scratch7/rob/dbg_index_tests/piscem/SRR1039508_mapped",
      "lib_type": "InwardUnstranded",
      "max_iter": 1500,
      "num_bootstraps": 0,
      "num_threads": 16,
      "output": "quant/SRR1039508"
    },
    "signatures": {
      "sha256_names": "e3f718f453cd3a749c2862a2c0a86ab6baf50529a5b73bf9cbfb95df31847542",
      "sha256_seqs": "fd88296600b98e1333273e657145662fd489a41f40296d387ab4265db2dc0f1c",
      "sha512_names": "df3e5b8133dc55ad4af4a9faab6bf34a031d3b8284b8f7e19526088a033cce1ddda1ebef5731fe4b094842ed48536b0095739af8b4ced 3fcfc897be0bdda5e23",
      "sha512_seqs": "373e06811ade1419f3640c476fbb8fd1e75b4eae2e43c189fb2f9da3b6f011d44cbf5027d9e48337a60fd2d6ec1d016fa9f23fce4bc01b f60690e81577d4c3de"
    }
  }


.. autosummary::
   :toctree: generated

    piscem-infer
