# SVELDT

SVELDT is bioinformatic tool that refinces SV using long read analysis

## Remarks

SVELDT is written in C++ and uses htslib library written in C. Please first install and compile htslib. Instructions are given below.

Please visit [htslib](https://github.com/samtools/htslib/tree/4604554d424406c6764af8da17b370c1b525ae1a) repository for detailed instructions.

## Installation

You need to first clone htslib and its submodule htscodecs. Please run following command inside the cloned repository of SVELDT.

```

    git submodule update --init --recursive

```

Then, go to the htslib directory to compile library.

```

    cd htslib

    ./configure
    make
    make install

```

The official installation documentation suggest to rum ./configure command if you need optional functionality to be enabled. If not, you can skip the command.

The default installation directory of libraries,library header files, utilities, several manual pages, and a pkgconfig is /usr/local.
This might cause permission errors if you installed the htslib in machines where you do not own the root privilege. Also, this configuration installs the library globally, but, you can change the directory by running following command 

```

    # instead of 
    make install
    
    # run     
    make prefix=DIR install

```

If after these, you should have installed and compiled htslib. If you encounter with any error, please visit [htslib](https://github.com/samtools/htslib/blob/4604554d424406c6764af8da17b370c1b525ae1a/INSTALL) and try to follow instructions there.
Now, you need to compile SVELDT. In order to do this please run following commands

```

    cd ..   # You need to go SVELDT directory
    
    make

```

# Usage

Great, now you can use the program by running 

```

    ./sveldt --bam bam_file --vcf vcf_file --output output_file_name(default_output.txt)

```

where bam_file, vcf_file and output_file_name are the input files.
