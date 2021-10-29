# dipdiff

This is a simple SV calling package for diploid assemblies. It uses a modified version of [svim-asm](https://github.com/eldariont/svim-asm).
The package includes its own version [minimap2](https://github.com/lh3/minimap2) to ensure reproducibility between runs, 
as the result might be dependent on the aligner version and parameters.

## Version 0.1

Installation & running
----------------------

```
git clone https://github.com/fenderglass/dipdiff
cd dipdiff
git submodule update --init --recursive
make
```

Afterwards, you can execute:

```
dipdiff.py
```

Output files
------------

The output directory will contain `variants.vcf` file with structural variants.


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


