import argparse
import gzip
import logging

import hail as hl

#from utils.clinvar import CLINVAR_GOLD_STARS_LOOKUP, download_and_import_latest_clinvar_vcf

from utils import (
    get_expr_for_alt_allele,
    get_expr_for_contig,
    get_expr_for_ref_allele,
    get_expr_for_start_pos,
    get_expr_for_variant_id,
    get_expr_for_xpos,
    get_expr_for_vep_consequence_terms_set,
    get_expr_for_vep_gene_id_to_consequence_map,
    get_expr_for_vep_gene_ids_set,
    get_expr_for_vep_protein_domains_set,
    get_expr_for_vep_sorted_transcript_consequences_array,
    get_expr_for_vep_transcript_id_to_consequence_map,
    get_expr_for_vep_transcript_ids_set,
    get_expr_for_worst_transcript_consequence_annotations_struct,
    get_expr_for_variant_ids,
)

from export_ht_to_es import *

logger = logging.getLogger()

CLINVAR_GOLD_STARS_LOOKUP = hl.dict(
    {
        "no_interpretation_for_the_single_variant": 0,
        "no_assertion_provided": 0,
        "no_assertion_criteria_provided": 0,
        "criteria_provided,_single_submitter": 1,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "reviewed_by_expert_panel": 3,
        "practice_guideline": 4,
    }
)



def _parse_clinvar_release_date(local_vcf_path: str) -> str:
    """Parse clinvar release date from the VCF header.

    Args:
        local_vcf_path (str): clinvar vcf path on the local file system.

    Returns:
        str: return VCF release date as string, or None if release date not found in header.
    """
    with gzip.open(local_vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##fileDate="):
                clinvar_release_date = line.split("=")[-1].strip()
                return clinvar_release_date

            if not line.startswith("#"):
                return None

    return None


'''
p = argparse.ArgumentParser()
p.add_argument("-g", "--genome-version", help="Genome build: 37 or 38", choices=["37", "38"], required=True)
p.add_argument("-H", "--host", help="Elasticsearch host or IP", required=True)
p.add_argument("-p", "--port", help="Elasticsearch port", default=9200, type=int)
p.add_argument("-i", "--index-name", help="Elasticsearch index name")
p.add_argument("-t", "--index-type", help="Elasticsearch index type", default="variant")
p.add_argument("-s", "--num-shards", help="Number of elasticsearch shards", default=1, type=int)
p.add_argument("-b", "--es-block-size", help="Elasticsearch block size to use when exporting", default=200, type=int)
args = p.parse_args()


if args.index_name:
    index_name = args.index_name.lower()
else:
    index_name = "clinvar_grch{}".format(args.genome_version)
'''

def import_vcf(
        vcf_path: str,
        genome_version: str,
        min_partitions: int = None,
        force_bgz: bool = True,
        drop_samples: bool = False,
        skip_invalid_loci: bool = False,
        split_multi_alleles: bool = True):
    """Import vcf and return MatrixTable.

    :param str vcf_path: MT to annotate with VEP
    :param str genome_version: "37" or "38"
    :param int min_partitions: min partitions
    :param bool force_bgz: read .gz as a bgzipped file
    :param bool drop_samples: if True, discard genotype info
    :param bool skip_invalid_loci: if True, skip loci that are not consistent with the reference_genome.
    """

    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {genome_version}")

    logger.info(f"\n==> import vcf: {vcf_path}")

    # add (or remove) "chr" prefix from vcf chroms so they match the reference
    ref = hl.get_reference(f"GRCh{genome_version}")
    contig_recoding = {
        **{ref_contig.replace("chr", ""): ref_contig for ref_contig in ref.contigs if "chr" in ref_contig},
        **{f"chr{ref_contig}": ref_contig for ref_contig in ref.contigs if "chr" not in ref_contig}}

    mt = hl.import_vcf(
        vcf_path,
        reference_genome=f"GRCh{genome_version}",
        contig_recoding=contig_recoding,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        drop_samples=drop_samples,
        skip_invalid_loci=skip_invalid_loci)

    mt = mt.annotate_globals(sourceFilePath=vcf_path, genomeVersion=genome_version)

    mt = mt.annotate_rows(
        original_alt_alleles=hl.or_missing(hl.len(mt.alleles) > 2, get_expr_for_variant_ids(mt.locus, mt.alleles))
    )

    if split_multi_alleles:
        mt = hl.split_multi_hts(mt)
        mt = mt.key_rows_by(**hl.min_rep(mt.locus, mt.alleles))

    return mt



def populate_clinvar():

    #clinvar_release_date = _parse_clinvar_release_date('clinvar.vcf.gz')
    #mt = import_vcf('clinvar.vcf.gz', "37", drop_samples=True, min_partitions=2000, skip_invalid_loci=True)
    #mt = mt.annotate_globals(version=clinvar_release_date)


    '''
    print("\n=== Running VEP ===")
    mt = hl.vep(mt, 'vep85-loftee-local.json', name="vep")

    print("\n=== Processing ===")
    mt = mt.annotate_rows(
        sortedTranscriptConsequences=get_expr_for_vep_sorted_transcript_consequences_array(vep_root=mt.vep)
    )

    mt = mt.annotate_rows(
        main_transcript=get_expr_for_worst_transcript_consequence_annotations_struct(
            vep_sorted_transcript_consequences_root=mt.sortedTranscriptConsequences
        )
    )

    mt = mt.annotate_rows(
        gene_ids=get_expr_for_vep_gene_ids_set(
            vep_transcript_consequences_root=mt.sortedTranscriptConsequences
        ),
    )

    review_status_str = hl.delimit(hl.sorted(hl.array(hl.set(mt.info.CLNREVSTAT)), key=lambda s: s.replace("^_", "z")))

    mt = mt.select_rows(
        allele_id=mt.info.ALLELEID,
        alt=get_expr_for_alt_allele(mt),
        chrom=get_expr_for_contig(mt.locus),
        clinical_significance=hl.delimit(hl.sorted(hl.array(hl.set(mt.info.CLNSIG)), key=lambda s: s.replace("^_", "z"))),
        domains=get_expr_for_vep_protein_domains_set(vep_transcript_consequences_root=mt.vep.transcript_consequences),
        gene_ids=mt.gene_ids,
        gene_id_to_consequence_json=get_expr_for_vep_gene_id_to_consequence_map(
            vep_sorted_transcript_consequences_root=mt.sortedTranscriptConsequences,
            gene_ids=mt.gene_ids
        ),
        gold_stars=CLINVAR_GOLD_STARS_LOOKUP[review_status_str],
        **{f"main_transcript_{field}": mt.main_transcript[field] for field in mt.main_transcript.dtype.fields},
        pos=get_expr_for_start_pos(mt),
        ref=get_expr_for_ref_allele(mt),
        review_status=review_status_str,
        transcript_consequence_terms=get_expr_for_vep_consequence_terms_set(
            vep_transcript_consequences_root=mt.sortedTranscriptConsequences
        ),
        transcript_ids=get_expr_for_vep_transcript_ids_set(
            vep_transcript_consequences_root=mt.sortedTranscriptConsequences
        ),
        transcript_id_to_consequence_json=get_expr_for_vep_transcript_id_to_consequence_map(
            vep_transcript_consequences_root=mt.sortedTranscriptConsequences
        ),
        variant_id=get_expr_for_variant_id(mt),
        xpos=get_expr_for_xpos(mt.locus),
    )

    #print("\n=== Summary ===")
    #hl.summarize_variants(mt)


    # Drop key columns for export
    rows = mt.rows()
    rows = rows.order_by(rows.variant_id).drop("locus", "alleles")
    rows.write('clinvar.ht',overwrite=True)
    '''
    print("\n=== Exporting to Elasticsearch ===")
    rows = hl.read_table('clinvar.ht')
    export_ht_to_es(rows, index_name = 'clinvar_grch37',index_type = 'variant')



if __name__ == "__main__":
    hl.init()
    populate_clinvar()
