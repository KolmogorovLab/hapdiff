# dipdiff

This is a simple SV calling package for diploid assemblies. It uses a modified version of [svim-asm](https://github.com/eldariont/svim-asm).
The package includes its own version [minimap2](https://github.com/lh3/minimap2) to ensure reproducibility between runs, 
as the result might be dependent on the aligner version and parameters.

## Version 0.2

Quick start
-----------


Dipdiff takes as input reference genome and a pair of haplotypes, and outputs
structural vaiant calls in VCF format. A recommended way to run is the Docker distribution.

Next steps assume that your `ref.fasta`, `hap_1.fasta` and `hap_2.fasta` are in the same directory,
which will also be used for dipdiff output. If it is not the case, you might need to bind additional 
directories using the Docker's `-v / --volume` argument. The number of threads (`-t` argument)
should be adjusted according to the available resources.


```
cd directory_with_input
DD_DIR=`pwd`
docker run -v $DD_DIR:$DD_DIR -u `id -u`:`id -g` mkolmogo/dipdiff:0.2 \
  dipdiff.py --reference $DD_DIR/ref.fasta --pat $DD_DIR/hap_1.fasta --mat $DD_DIR/hap_2.fasta --out-dir $DD_DIR/dipdiff -t 20
```

Output files
------------

The output directory will contain `variants.vcf` file with structural variants.


Source Installation
-------------------

Alernatively, you can run dipdiff locally as follows.

```
git clone https://github.com/fenderglass/dipdiff
cd dipdiff
git submodule update --init
make
pip install -r requirements.txt
```

In addition, dipdiff requires [samtools](https://github.com/samtools) to be installed in your system.

Afterwards, you can execute:

```
./dipdiff.py --reference ref.fasta --pat hap_1.fasta --mat hap_2.fasta --out-dir out_path -t 20
```

Acknowledgements
----------------

The major parts of the dipdiff pipeline are:

* [minimap2](https://github.com/lh3/minimap2)
* [svim-asm](https://github.com/eldariont/svim-asm)


Authors
-------

The pipeline was developed at [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/), Benedict Paten's lab.

Main code contributors:
* Mikhail Kolmogorov


License
-------

dipdiff is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.
Other software included in this discrubution is released under either MIT or BSD licenses.


How to get help
---------------
A preferred way report any problems or ask questions is the 
[issue tracker](https://github.com/fenderglass/dipdiff/issues). 


