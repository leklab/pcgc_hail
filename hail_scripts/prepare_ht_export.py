import hail as hl

import pprint
import copy
import itertools
from collections import defaultdict, namedtuple, OrderedDict
from typing import *


GROUPS = ['adj', 'raw']
SEXES = ['male', 'female']
PROBAND = ['proband']
#POPS = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas']
POPS = ['afr', 'amr', 'eas', 'eur', 'oth', 'sas']

SORT_ORDER = ['popmax', 'group', 'pop', 'proband', 'subpop', 'sex']


pop_names = {
    'oth': 'Other',
    'afr': 'African-American/African',
    'amr': 'Latino',
    'eas': 'East Asian',
    'fin': 'Finnish',
    'eur': 'European',
    'nfe': 'Non-Finnish European',
    'sas': 'South Asian',
    'mde': 'Middle Eastern',
    'asj': 'Ashkenazi Jewish',
    'uniform': 'Uniform',
    'sas_non_consang': 'South Asian (F < 0.05)',
    'consanguineous': 'South Asian (F > 0.05)',
    'exac': 'ExAC',
    'bgr': 'Bulgarian (Eastern European)',
    'deu': 'German',
    'est': 'Estonian',
    'esp': 'Spanish',
    'gbr': 'British',
    'nwe': 'North-Western European',
    'seu': 'Southern European',
    'ita': 'Italian',
    'swe': 'Swedish',
    'chn': 'Chinese',
    'kor': 'Korean',
    'hkg': 'Hong Kong',
    'sgp': 'Singaporean',
    'twn': 'Taiwanese',
    'jpn': 'Japanese',
    'oea': 'Other East Asian',
    'oeu': 'Other European',
    'onf': 'Other Non-Finnish European',
    'unk': 'Unknown'
}

INFO_DICT = {
    'FS': {"Description": "Phred-scaled p-value of Fisher's exact test for strand bias"},
    'InbreedingCoeff': {
        "Description": "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"},
    'MQ': {"Description": "Root mean square of the mapping quality of reads across all samples"},
    'MQRankSum': {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"},
    'QD': {"Description": "Variant call confidence normalized by depth of sample reads supporting a variant"},
    'ReadPosRankSum': {"Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias"},
    'SOR': {"Description": "Strand bias estimated by the symmetric odds ratio test"},
    'VQSR_POSITIVE_TRAIN_SITE': {"Description": "Variant was used to build the positive training set of high-quality variants for VQSR"},
    'VQSR_NEGATIVE_TRAIN_SITE': {
        "Description": "Variant was used to build the negative training set of low-quality variants for VQSR"},
    'BaseQRankSum': {"Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities"},
    'ClippingRankSum': {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference number of hard clipped bases"},
    'DP': {"Description": "Depth of informative coverage for each sample; reads with MQ=255 or with bad mates are filtered"},
    'VQSLOD': {
        "Description": "Log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model"},
    'VQSR_culprit': {"Description": "Worst-performing annotation in the VQSR Gaussian mixture model"},
}






def make_label_combos(label_groups: Dict[str, List[str]]) -> List[str]:
    '''
    Make combinations of all possible labels for a supplied dictionary of label groups
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :return: list of all possible combinations of values for the supplied label groupings
    :rtype: list[str]
    '''
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: SORT_ORDER.index(x))[0]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x,y in itertools.product(anchor_val, make_label_combos(copy_label_groups)):
        combos.append('{0}_{1}'.format(x,y))
    return combos


def make_combo_header_text(preposition, group_types, combo_fields, prefix, faf=False):
    '''
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset
    :param str preposition: Relevant preposition to precede automatically generated text
    :param list of str group_types: List of grouping types, e.g. "sex" or "pop"
    :param list of str combo_fields: List of the specific values for each grouping type, for which the text is being generated
    :param str prefix: Subset of gnomAD
    :param bool faf: If True, use alternate logic to automatically populate descriptions for filter allele frequency annotations
    :return: String with automatically generated description text for a given set of combo fields
    :rtype: str
    '''
    combo_dict = dict(zip(group_types, combo_fields))

    if not faf:
        header_text = " " + preposition
        if 'sex' in combo_dict.keys():
            header_text = header_text + " " + combo_dict['sex']
        header_text = header_text + " samples"
        if 'subpop' in combo_dict.keys():
            header_text = header_text + f" of {pop_names[combo_dict['subpop']]} ancestry"
            combo_dict.pop('pop')
        if 'pop' in combo_dict.keys():
            header_text = header_text + f" of {pop_names[combo_dict['pop']]} ancestry"
        if prefix != 'gnomad':
            header_text = header_text + f" in the {prefix} subset"
        if 'group' in group_types:
            if combo_dict['group'] == 'raw':
                header_text = header_text + ", before removing low-confidence genotypes"
    else:
        header_text = ""
        if 'pop' in combo_dict.keys():
            header_text = f" of {pop_names[combo_dict['pop']]} ancestry"
        if prefix != 'gnomad':
            header_text = header_text + f" in the {prefix} subset"
    return header_text


def make_info_dict(prefix, label_groups=None, bin_edges=None, faf=False, popmax=False, age_hist_data=None):
    '''
    Generate dictionary of Number and Description attributes to be used in the VCF header
    :param str prefix: Subset of gnomAD
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param dict bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation
    :param bool faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations
    :param bool popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations
    :param str age_hist_data: Pipe-delimited string of age histograms, from get_age_distributions (somewhat required if popmax == True)
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    '''
    info_dict = dict()

    if popmax:
        popmax_text = "" if prefix == 'gnomad' else f" in the {prefix} subset"
        popmax_dict = {
            f'{prefix}_popmax': {"Number": "A",
                                 "Description": "Population with maximum AF{}".format(popmax_text)},
            f'{prefix}_AC_popmax': {"Number": "A",
                                    "Description": "Allele count in the population with the maximum AF{}".format(popmax_text)},
            f'{prefix}_AN_popmax': {"Number": "A",
                                    "Description": "Total number of alleles in the population with the maximum AF{}".format(popmax_text)},
            f'{prefix}_AF_popmax': {"Number": "A",
                                    "Description": "Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry){}".format(popmax_text)},
            f'{prefix}_nhomalt_popmax': {"Number": "A",
                                         "Description": "Count of homozygous individuals in the population with the maximum allele frequency{}".format(popmax_text)}
        }
        info_dict.update(popmax_dict)
        if prefix == 'gnomad':
            age_hist_dict = {
                f"{prefix}_age_hist_het_bin_freq": {"Number": "A",
                                                    "Description": f"Histogram of ages of heterozygous individuals; bin edges are: {bin_edges[f'{prefix}_het']}; total number of individuals of any genotype bin: {age_hist_data}"},
                f"{prefix}_age_hist_het_n_smaller": {"Number": "A",
                                                     "Description": "Count of age values falling below lowest histogram bin edge for heterozygous individuals"},
                f"{prefix}_age_hist_het_n_larger": {"Number": "A",
                                                    "Description": "Count of age values falling above highest histogram bin edge for heterozygous individuals"},
                f"{prefix}_age_hist_hom_bin_freq": {"Number": "A",
                                                    "Description": f"Histogram of ages of homozygous alternate individuals; bin edges are: {bin_edges[f'{prefix}_hom']}; total number of individuals of any genotype bin: {age_hist_data}"},
                f"{prefix}_age_hist_hom_n_smaller": {"Number": "A",
                                                     "Description": "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals"},
                f"{prefix}_age_hist_hom_n_larger": {"Number": "A",
                                                    "Description": "Count of age values falling above highest histogram bin edge for homozygous alternate individuals"}
            }
            info_dict.update(age_hist_dict)
    else:
        group_types = sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
        combos = make_label_combos(label_groups)

        for combo in combos:
            combo_fields = combo.split("_")
            if not faf:
                combo_dict = {
                    f"{prefix}_AC_{combo}": {"Number": "A",
                                             "Description": "Alternate allele count{}".format(make_combo_header_text('for', group_types, combo_fields, prefix))},
                    f"{prefix}_AN_{combo}": {"Number": "1",
                                             "Description": "Total number of alleles{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))},
                    f"{prefix}_AF_{combo}": {"Number": "A",
                                             "Description": "Alternate allele frequency{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))},
                    f"{prefix}_nhomalt_{combo}": {"Number": "A",
                                                  "Description": "Count of homozygous individuals{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))}
                }
            else:
                combo_dict = {
                    f"{prefix}_faf95_{combo}": {"Number": "A",
                                                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples{}".format(make_combo_header_text(None, group_types, combo_fields, prefix, faf=True))},
                    f"{prefix}_faf99_{combo}": {"Number": "A",
                                                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples{}".format(make_combo_header_text(None, group_types, combo_fields, prefix, faf=True))}
                }
            info_dict.update(combo_dict)
    return info_dict



def make_info_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    '''
    Make Hail expression for variant annotations to be included in VCF INFO
    :param Table ht: Table containing variant annotations to be reformatted for VCF export
    :return: Dictionary containing Hail expressions for relevant INFO annotations
    :rtype: Dict of str: Expression
    '''
    info_dict = {
        'FS': ht.info.FS,
        'InbreedingCoeff': ht.info.InbreedingCoeff,
        'MQ': ht.info.MQ,
        'MQRankSum': ht.info.MQRankSum,
        'QD': ht.info.QD,
        'ReadPosRankSum': ht.info.ReadPosRankSum,
        'SOR': ht.info.SOR,
        'VQSR_POSITIVE_TRAIN_SITE': ht.info.POSITIVE_TRAIN_SITE,
        'VQSR_NEGATIVE_TRAIN_SITE': ht.info.NEGATIVE_TRAIN_SITE,
#        'BaseQRankSum': ht.allele_info.BaseQRankSum,
#        'ClippingRankSum': ht.allele_info.ClippingRankSum,
#        'DP': ht.allele_info.DP,
#        'VQSLOD': ht.allele_info.VQSLOD,
#        'VQSR_culprit': ht.allele_info.culprit,
    }

    return info_dict

def index_globals(freq_meta, label_groups):
    '''
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry
    in the freq_meta array annotation
    :param list of str freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :return: Dictionary keyed by specified label grouping combinations, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    combos = make_label_combos(label_groups)
    index_dict = {}

    for combo in combos:
        combo_fields = combo.split("_")
        for i,v in enumerate(freq_meta):
            if set(v.values()) == set(combo_fields):
                index_dict.update({f"{combo}": i})
    return index_dict


def make_freq_meta_index_dict(freq_meta):
    '''
    Make dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices
    :param list of str freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    #pprint.pprint(freq_meta)
    index_dict = index_globals(freq_meta, dict(group=GROUPS))
    #pprint.pprint(index_dict)
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS)))
    #pprint.pprint(index_dict)
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, proband=PROBAND)))
    #pprint.pprint(index_dict)

    #index_dict.update(index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)))
    #index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)))
    #index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['nfe'], subpop=NFE_SUBPOPS)))
    #index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['eas'], subpop=EAS_SUBPOPS)))
    return index_dict

def make_index_dict(ht):
    '''
    Create a look-up Dictionary for entries contained in the frequency annotation array
    :param Table ht: Table containing freq_meta global annotation to be indexed
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    freq_meta = hl.eval(ht.globals.freq_meta)
    index_dict = make_freq_meta_index_dict(freq_meta)
    return index_dict

def unfurl_nested_annotations(ht):
    '''
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values
    :param Table ht: Hail Table containing the nested variant annotation arrays to be unfurled
    :return: Dictionary containing variant annotations and their corresponding values
    :rtype: Dict of str: Expression
    '''
    expr_dict = dict()
    freq_index_dict = make_index_dict(ht)

    #ML: Removing adj prefix so only shows population and adj, raw
    #pprint.pprint(freq_index_dict)
    new_freq_index_dict = {i.replace('adj_', 'gnomad_'): j for i,j in freq_index_dict.items()}
    #pprint.pprint(new_freq_index_dict)

    #for k, i in hl.eval(ht.globals.freq_index_dict).items():
    for k, i in new_freq_index_dict.items():
        entry = k.split("_")
        if entry[0] == "non":
            prefix = "_".join(entry[:2])
            combo_fields = ['adj'] + entry[2:]
            if combo_fields == ['adj', 'raw']:
                combo_fields = ['raw']
        
        elif entry[0] == "raw":
            combo_fields = ['raw']

        else:
            prefix = entry[0]
            combo_fields = ['adj'] + entry[1:]
            if combo_fields == ['adj', 'raw']:
                combo_fields = ['raw']

        combo = "_".join(combo_fields)

        combo_dict = {
            f"AC_{combo}": ht.freq[i].AC,
            f"AN_{combo}": ht.freq[i].AN,
            f"AF_{combo}": ht.freq[i].AF,
            f"nhomalt_{combo}": ht.freq[i].homozygote_count
        }

        #combo_dict = {
        #    f"{prefix}_AC_{combo}": ht.freq[i].AC,
        #    f"{prefix}_AN_{combo}": ht.freq[i].AN,
        #    f"{prefix}_AF_{combo}": ht.freq[i].AF,
        #    f"{prefix}_nhomalt_{combo}": ht.freq[i].homozygote_count
        #}

        expr_dict.update(combo_dict)

    #pprint.pprint(expr_dict)
    return expr_dict



def prepare_ht_export(ht: hl.Table) -> hl.Table:

    subset_list = ['gnomad']

    for subset in subset_list:
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS)))
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, pop=POPS)))

    new_info_dict = {i.replace('gnomad_', '').replace('_adj', ''): j for i,j in INFO_DICT.items()}

    ht = ht.annotate(info=hl.struct(**make_info_expr(ht)))
    ht = ht.annotate(info=ht.info.annotate(**unfurl_nested_annotations(ht)))


    #ht = ht.select('info', 'filters', 'rsid', 'qual','vep')
    ht = ht.select('info', 'filters', 'rsid', 'qual')


    header_dict = {'info': new_info_dict}
                   #'filter': make_filter_dict(ht)}

    return ht



