import pprint
from annotate_frequencies import *
from generate_split_alleles import *
from prepare_ht_export import *


def run_pipeline():
    hl.init(log='./hail_annotation_pipeline.log')
    mt = hl.import_vcf('vcf_files/pcgc_chr20_slice.vcf.bgz',reference_genome='GRCh37')

    #Split alleles
    mt = generate_split_alleles(mt)
    #pprint.pprint(mt.describe())
    #pprint.pprint(mt.show(include_row_fields=True))

    #Annotate Population frequencies for now
    meta_ht = hl.import_table('vcf_files/pcgc_meta.tsv',delimiter='\t',key='ID')
    ht = annotate_frequencies(mt,meta_ht)
    #pprint.pprint(ht.describe())
    #pprint.pprint(ht.show())

    #VEP Annotate the Hail table (ie. sites-only)
    ht = hl.vep(ht, 'vcf_files/vep85-loftee-local.json')
    #pprint.pprint(ht.describe())
    #pprint.pprint(ht.show())

    ht = prepare_ht_export(ht)
    pprint.pprint(ht.describe())
    pprint.pprint(ht.show())


if __name__ == '__main__':
    run_pipeline()

