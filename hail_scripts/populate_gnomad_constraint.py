import argparse
import pprint
import hail as hl

from export_ht_to_es import *

#gsutil -m cp -r gs://gnomad-public/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.lof_metrics.by_transcript.ht .
#gsutil cp gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/constraint_final_standard.txt.bgz .

def populate_constraint():

    #ds = hl.read_table('gnomad.v2.1.1.lof_metrics.by_transcript.ht')
    #ds = hl.import_table('constraint_final_standard.txt.bgz',delimiter='\t',key='transcript',impute=True)
    ds = hl.import_table('constraint_final_cleaned.txt.bgz',delimiter='\t',key='transcript',impute=True)

    #ds = hl.import_table('missing_small.txt',delimiter='\t',key='transcript',impute=True)

    # The globals in the Hail table cause a serialization error during Elasticsearch export
    ds = ds.select_globals()
    pprint.pprint(ds.describe())
    '''
    population_dict_fields = [
        "pop_no_lofs",
        "pop_obs_het_lof",
        "pop_obs_hom_lof",
        "pop_defined",
        "pop_p",
    ]

    populations = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]

    # Convert dicts to structs for Elasticsearch export
    ds = ds.annotate(
        **{
            f: hl.struct(**{pop: ds[f][pop] for pop in populations})
            for f in population_dict_fields
        }
    )

    '''

    # Convert interval to struct for Elasticsearch export
    '''
    ds = ds.annotate(
        interval=hl.struct(
            chrom=ds.interval.start.contig,
            start=ds.interval.start.position,
            end=ds.interval.end.position,
        )
    )
    
    ds = ds.key_by()
    '''
    ds = ds.transmute(gene_name=ds.gene, transcript_id=ds.transcript)

    #ds.write(args.output_url)

    '''
    ds = ds.select('exp_lof','exp_mis','exp_syn','obs_lof','obs_mis','obs_syn',
                    'oe_lof','oe_lof_lower','oe_mis','oe_mis_lower','oe_mis_upper',
                    'oe_syn','oe_syn_lower','oe_syn_upper',
                    'lof_z','mis_z','syn_z',
                    'pLI','pNull','pRec')
    '''
    '''
    ds = ds.select('exp_lof','exp_mis','exp_syn','obs_lof','obs_mis','obs_syn',
                    'lof_z','mis_z','syn_z',
                    'pLI','pNull','pRec')

    '''
    pprint.pprint(ds.describe())
    pprint.pprint(ds.show())

    export_ht_to_es(ds, index_name = 'gnomad_constraint_2_1_1',index_type = 'constraint')



if __name__ == "__main__":
    hl.init()
    populate_constraint()
