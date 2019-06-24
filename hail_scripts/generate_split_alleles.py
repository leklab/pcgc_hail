import hail as hl

def add_variant_type(alt_alleles: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """
    Get Struct of variant_type and n_alt_alleles from ArrayExpression of Strings (all alleles)
    """
    ref = alt_alleles[0]
    alts = alt_alleles[1:]
    non_star_alleles = hl.filter(lambda a: a != '*', alts)
    return hl.struct(variant_type=hl.cond(
        hl.all(lambda a: hl.is_snp(ref, a), non_star_alleles),
        hl.cond(hl.len(non_star_alleles) > 1, "multi-snv", "snv"),
        hl.cond(
            hl.all(lambda a: hl.is_indel(ref, a), non_star_alleles),
            hl.cond(hl.len(non_star_alleles) > 1, "multi-indel", "indel"),
            "mixed")
    ), n_alt_alleles=hl.len(non_star_alleles))


def generate_split_alleles(mt: hl.MatrixTable) -> hl.Table:

    allele_data = hl.struct(nonsplit_alleles=mt.alleles,
                            has_star=hl.any(lambda a: a == '*', mt.alleles))

    mt = mt.annotate_rows(allele_data=allele_data.annotate(**add_variant_type(mt.alleles)))
    mt = hl.split_multi_hts(mt,left_aligned=True)

    allele_type = (hl.case()
                   .when(hl.is_snp(mt.alleles[0], mt.alleles[1]), 'snv')
                   .when(hl.is_insertion(mt.alleles[0], mt.alleles[1]), 'ins')
                   .when(hl.is_deletion(mt.alleles[0], mt.alleles[1]), 'del')
                   .default('complex')
                   )
    mt = mt.annotate_rows(allele_data=mt.allele_data.annotate(allele_type=allele_type,
                                                              was_mixed=mt.allele_data.variant_type == 'mixed'))
    return mt

