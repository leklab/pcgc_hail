import hail as hl
import pprint
from export_ht_to_es import *

tissue_abbr = {
	'Adipose - Subcutaneous' : 'adiposeSubcutaneous',
	'Adipose - Visceral (Omentum)' : 'adiposeVisceralOmentum',
	'Adrenal Gland' : 'adrenalGland',
	'Artery - Aorta' : 'arteryAorta',
	'Artery - Coronary' : 'arteryCoronary',
	'Artery - Tibial' : 'arteryTibial',
	'Bladder' : 'bladder',
	'Brain - Amygdala' : 'brainAmygdala',
	'Brain - Anterior cingulate cortex (BA24)' : 'brainAnteriorcingulatecortexBa24',
	'Brain - Caudate (basal ganglia)' : 'brainCaudateBasalganglia',
	'Brain - Cerebellar Hemisphere' : 'brainCerebellarhemisphere',
	'Brain - Cerebellum' : 'brainCerebellum',
	'Brain - Cortex' : 'brainCortex',
	'Brain - Frontal Cortex (BA9)' : 'brainFrontalcortexBa9',
	'Brain - Hippocampus' : 'brainHippocampus',
	'Brain - Hypothalamus' : 'brainHypothalamus',
	'Brain - Nucleus accumbens (basal ganglia)' : 'brainNucleusaccumbensBasalganglia',
	'Brain - Putamen (basal ganglia)' : 'brainPutamenBasalganglia',
	'Brain - Spinal cord (cervical c-1)' : 'brainSpinalcordCervicalc1',
	'Brain - Substantia nigra' : 'brainSubstantianigra',
	'Breast - Mammary Tissue' : 'breastMammarytissue',
	'Cells - EBV-transformed lymphocytes' : 'cellsEbvTransformedlymphocytes',
	'Cells - Transformed fibroblasts' : 'cellsTransformedfibroblasts',
	'Cervix - Ectocervix' : 'cervixEctocervix',
	'Cervix - Endocervix' : 'cervixEndocervix',
	'Colon - Sigmoid' : 'colonSigmoid',
	'Colon - Transverse' : 'colonTransverse',
	'Esophagus - Gastroesophageal Junction' : 'esophagusGastroesophagealjunction',
	'Esophagus - Mucosa' : 'esophagusMucosa',
	'Esophagus - Muscularis' : 'esophagusMuscularis',
	'Fallopian Tube' : 'fallopianTube',
	'Heart - Atrial Appendage' : 'heartAtrialappendage',
	'Heart - Left Ventricle' : 'heartLeftventricle',
	'Kidney - Cortex' : 'kidneyCortex',
	'Liver' : 'liver',
	'Lung' : 'lung',
	'Minor Salivary Gland' : 'minorSalivaryGland',
	'Muscle - Skeletal' : 'muscleSkeletal',
	'Nerve - Tibial' : 'nerveTibial',
	'Ovary' : 'ovary',
	'Pancreas' : 'pancreas',
	'Pituitary' : 'pituitary',
	'Prostate' : 'prostate',
	'Skin - Not Sun Exposed (Suprapubic)' : 'skinNotsunexposedSuprapubic',
	'Skin - Sun Exposed (Lower leg)' : 'skinSunexposedLowerleg',
	'Small Intestine - Terminal Ileum' : 'smallIntestineTerminalileum',
	'Spleen' : 'spleen',
	'Stomach' : 'stomach',
	'Testis' : 'testis',
	'Thyroid' : 'thyroid',
	'Uterus' : 'uterus',
	'Vagina' : 'vagina',
	'Whole Blood' : 'wholeBlood'
}



def populate_gtex():
	meta_ht = hl.import_table('/home/ml2529/gtex_data/GTEx_v7_Annotations_SampleAttributesDS.txt',delimiter='\t',key='SAMPID')
	mt = hl.import_matrix_table('/home/ml2529/gtex_data/ENSG00000177732.tsv', row_key='transcript_id', row_fields={'transcript_id': hl.tstr, 'gene_id': hl.tstr},entry_type=hl.tfloat32)
	#mt = hl.import_matrix_table('/home/ml2529/gtex_data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.bgz', row_key='transcript_id', row_fields={'transcript_id': hl.tstr, 'gene_id': hl.tstr},entry_type=hl.tfloat32)

	mt = mt.rename({'transcript_id': 'transcriptId', 'gene_id': 'geneId'})

	#pprint.pprint(meta_ht.describe())
	#pprint.pprint(gtex_mt.describe())

	mt = mt.annotate_cols(tissue = meta_ht[mt.col_id].SMTSD)

	#pprint.pprint(mt.describe())
	#pprint.pprint(mt.show(include_row_fields=True))

	cut_dict = {'tissue': hl.agg.filter(hl.is_defined(mt.tissue), hl.agg.counter(mt.tissue))}
	#pprint.pprint(cut_dict)

	cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
	#pprint.pprint(cut_data.tissue)

	#call_stats = hl.agg.filter(mt.tissue == 'Lung', hl.agg.mean(mt.x))
	#pprint.pprint(call_stats)

	#mt = mt.annotate_rows(Lung=call_stats)
	#pprint.pprint(mt.show(include_row_fields=True))

	for x in sorted(cut_data['tissue'].keys()):
		#pprint.pprint(x)
		call_stats = hl.agg.filter(mt.tissue == x, hl.agg.mean(mt.x))
		mt = mt.transmute_rows(**{f"{tissue_abbr[x]}": call_stats})

	#pprint.pprint(mt.show(include_row_fields=True))	

	ht = mt.rows()

	#ht.write('gtex_expression.ht',overwrite=True)

	export_ht_to_es(ht, index_name = 'gtex_tissue_tpms_by_transcript',index_type = 'tissue_tpms')

	'''
	sample_group_filters = [({}, True)]
	sample_group_filters.extend([({'tissue': tissue}, mt.tissue == tissue) for tissue in cut_data.tissue])
	mt = mt.select_cols(group_membership=tuple(x[1] for x in sample_group_filters))

	tissue_expression = []
	meta_expressions = []

	for i in range(len(sample_group_filters)):
		subgroup_dict = sample_group_filters[i][0]
		subgroup_dict['group'] = 'adj'
		#call_stats = hl.agg.filter(mt.group_membership[i], hl.agg.stats(mt.x))
		call_stats = hl.agg.filter(mt.group_membership[i], hl.agg.mean(mt.x))
		#pprint.pprint(subgroup_dict)
		pprint.pprint(call_stats)

		#call_stats_bind = hl.bind(lambda cs: cs.annotate(mean=cs.mean), call_stats)
		#tissue_expression.append(call_stats_bind)
		tissue_expression.append(call_stats)
		meta_expressions.append(subgroup_dict)



	global_expression = {
        'tissue_meta': meta_expressions
    }

	mt = mt.annotate_rows(tissue_expression=tissue_expression)	
	mt = mt.annotate_globals(**global_expression)


	mt = mt.rows()

	pprint.pprint(mt.describe())
	#pprint.pprint(mt.show(include_row_fields=True))
	pprint.pprint(mt.show())
	'''


if __name__ == "__main__":
	hl.init()
	populate_gtex()
