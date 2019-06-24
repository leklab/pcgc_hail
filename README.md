# pcgc_hail
Hail v0.2 pipeline used to VCF file from PCGC for export into elastic search.

## Requirements
Python 3.6
Spark-2.4

## Example
python submit.py --run-locally hail_annotate_pipeline.py --spark-home $SPARK_HOME --driver-memory 16G --executor-memory 8G
