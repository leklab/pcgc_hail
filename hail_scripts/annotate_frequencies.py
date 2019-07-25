import hail as hl
from typing import *
import pprint

def get_adj_expr(
        gt_expr: hl.expr.CallExpression,
        gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        ad_expr: hl.expr.ArrayNumericExpression,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.expr.BooleanExpression:
    """
    Gets adj genotype annotation.
    Defaults correspond to gnomAD values.
    """
    return (
            (gq_expr >= adj_gq) &
            hl.cond(
                gt_expr.is_haploid(),
                dp_expr >= haploid_adj_dp,
                dp_expr >= adj_dp
            ) &
            (
                hl.case()
                .when(~gt_expr.is_het(), True)
                .when(gt_expr.is_het_ref(), ad_expr[1] / dp_expr >= adj_ab)
                .default((ad_expr[0] / dp_expr >= adj_ab ) & (ad_expr[1] / dp_expr >= adj_ab ))
            )
    )

def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp))


def annotate_frequencies(mt: hl.MatrixTable, meta_ht: hl.Table) -> hl.Table:

    #meta_ht = hl.import_table('vcf_files/pcgc_meta.tsv',delimiter='\t',key='ID')
    #mt = hl.import_vcf('vcf_files/pcgc_chr20_slice.vcf.bgz',reference_genome='GRCh37')

    mt = mt.annotate_cols(pop = meta_ht[mt.s].Ethnicity)
    mt = mt.annotate_cols(proband = meta_ht[mt.s].Proband)


    mt = annotate_adj(mt)


    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.pop), hl.agg.counter(mt.pop))}
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend([
        ({'pop': pop}, mt.pop == pop) for pop in cut_data.pop])

    sample_group_filters.extend([({'proband': 'proband'}, mt.proband == 'Yes')])

    mt = mt.select_cols(group_membership=tuple(x[1] for x in sample_group_filters))

    frequency_expression = []
    meta_expressions = []

    for i in range(len(sample_group_filters)):
        subgroup_dict = sample_group_filters[i][0]
        #pprint.pprint(subgroup_dict)
        subgroup_dict['group'] = 'adj'

        call_stats = hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles))

        call_stats_bind = hl.bind(lambda cs: cs.annotate(
            AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
        ), call_stats)
        frequency_expression.append(call_stats_bind)
        meta_expressions.append(subgroup_dict)


    raw_stats = hl.agg.call_stats(mt.GT, mt.alleles)

    raw_stats_bind = hl.bind(lambda cs: cs.annotate(
        AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
    ), raw_stats)

    frequency_expression.insert(1, raw_stats_bind)
    meta_expressions.insert(1, {'group': 'raw'})

    '''
    proband_stats = hl.agg.filter(mt.proband == 'Yes' & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles))

    proband_stats_bind = hl.bind(lambda cs: cs.annotate(
        AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
    ), proband_stats)

    frequency_expression.insert(1, proband_stats_bind)
    meta_expressions.insert(1, {'group': 'proband'})
    '''


    print(f'Calculating {len(frequency_expression)} aggregators...')

    global_expression = {
        'freq_meta': meta_expressions
    }

    mt = mt.annotate_rows(freq=frequency_expression)
    mt = mt.annotate_globals(**global_expression)

    return mt.rows()




