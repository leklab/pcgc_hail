# pcgc_hail
Hail v0.2 pipeline used to VCF file from PCGC for export into elastic search.

## Requirements
```
Python 3.6
Java 1.8
Spark-2.4
VEP v85 (with LOFETEE plugin)
```

## High level summary of Pipeline steps
```
Step 1: Break multi-allelics to bi-allelics and use minimal representation of alleles
Step 2: Count of alleles by population
Step 3: Annotate using Variant Effect Predictor (VEP) with LOFTEE plugin
Step 4: Reformat for export to Hail Table and VCF file
Step 5: Reformat and export to Elastic Search
```

## Example Usage
```
python submit.py --run-locally ./hail_scripts/hail_annotate_pipeline.py --spark-home $SPARK_HOME --driver-memory 16G --executor-memory 8G -i input_file -m meta_file
```

## Wookie mistakes
Python 3.6 is not the default python

$SPARK_HOME is set to an older version of Spark

