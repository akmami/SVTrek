# SVTrek

SVTrek is a bioinformatics tool for evaluating structural variation (SV) discoveries and can also perform SV discovery independently using long-read analysis.

## Remarks

SVTrek is written in C and utilizes the `HTSlib` library.

## Installation

Follow these steps to install `SVTrek`:

```
# clone the repository
git clone --recursive --depth 1 https://github.com/akmami/SVTrek.git

# install submodules
cd SVTrek
make install

# compile the program
make
```

## Usage
```
./svtrek [-b|--bam BAM] [-v|--vcf VCF file] [OPTIONS]
```

## Required Parameters
- `-b, --bam <BAM>`
  - Specifies the BAM file to be processed.
- `-v, --vcf <VCF file>`
  - Specifies the VCF file to be used.

## Options
- `-o, --output <filename>`
  - Specifies the output filename.
  - **Default:** `svtrek.out`
- `-t <num>`
  - Number of threads to use for processing.
  - **Default:** `4`
- `--verbose`
  - Enables verbose output.
  - **Default:** `false`
- `--wider-interval <num>`
  - Defines the offset interval for the start of the reads (DEL-START).
  - **Default:** `20000`
- `--median-interval <num>`
  - Defines the offset interval for the start of the reads (INS).
  - **Default:** `10000`
- `--narrow-interval <num>`
  - Defines the offset interval for the end of the reads (DEL-END).
  - **Default:** `2000`
- `--consensus-interval-range <num>`
  - Specifies the offset that determines whether to consider the locations into the refinement.
  - **Default:** `500`
- `--consensus-interval <num>`
  - Specifies the interval that determines whether reads are considered to be in the same position.
  - **Default:** `5`
- `--consensus-min-count <num>`
  - Minimum number of elements required for consensus determination.
  - **Default:** `3`

### Example Usage
```
./svtrek -b input.bam -v input.vcf
```

