import pprint
import argparse

from annotate_frequencies import *
from generate_split_alleles import *
from prepare_ht_export import *
from prepare_ht_for_es import *
from export_ht_to_es import *


def run_pipeline(args):
    hl.init(log='./hail_annotation_pipeline.log')

    #mt = hl.import_vcf('vcf_files/pcgc_chr20_slice.vcf.bgz',reference_genome='GRCh37')
    mt = hl.import_vcf(args.vcf,reference_genome='GRCh37')

    #Split alleles
    mt = generate_split_alleles(mt)
    #pprint.pprint(mt.describe())
    #pprint.pprint(mt.show(include_row_fields=True))

    #Annotate Population frequencies for now
    meta_ht = hl.import_table(args.meta,delimiter='\t',key='ID')
    ht = annotate_frequencies(mt,meta_ht)
    #pprint.pprint(ht.describe())
    #pprint.pprint(ht.show())

    #VEP Annotate the Hail table (ie. sites-only)
    #ht = hl.vep(ht, 'vep85-loftee-local.json')
    #pprint.pprint(ht.describe())
    #pprint.pprint(ht.show())

    ht = prepare_ht_export(ht)
    #pprint.pprint(ht.describe()) 
    #pprint.pprint(ht.show())

    ht = prepare_ht_for_es(ht)
    #pprint.pprint(ht.describe())
    #pprint.pprint(ht.show())

    ht.write('pcgc_chr20_100samples.ht',overwrite=True)
    #export_ht_to_es(ht)

    #ht = hl.read_table('/home/ml2529/PCGC_dev/data/pcgc_chr20_100samples.ht')
    #ht = hl.read_table('/home/ml2529/PCGC_dev/data/pcgc_exomes.ht')
    #export_ht_to_es(ht)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='bgzipped VCF file (.vcf.bgz)', required=True)
    parser.add_argument('--meta', '-m', help='Meta file containing sample population and sex', required=True)

    args = parser.parse_args()
    run_pipeline(args)
