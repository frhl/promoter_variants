#!/usr/bin/env python3

import hail as hl
import os
import argparse

def hail_init(chrom=None, log_prefix='get_vcf'):
    r'''Initialize Hail '''
    assert chrom in range(1,23), 'only autosomes allowed'
    n_slots = os.environ.get('NSLOTS', 1)
    chr_suffix = '' if chr is None else f'_chr{chrom}'
    WD = '/well/lindgren/UKBIOBANK/flassen/projects/promoter_variants'
    hl.init(log=f'{WD}/logs/hail-{log_prefix}{chr_suffix}.log',
            default_reference='GRCh37',
            master=f'local[{n_slots}]')

def get_table(input_path, input_type, cache=False):
    r'''Import mt/vcf/plink tables '''
    if input_type=='mt':
        mt = hl.read_matrix_table(input_path)
    elif input_type=='vcf':
        mt = hl.import_vcf(input_path, force_bgz=True, array_elements_required=False)
    elif input_type=='plink':
        mt = hl.import_plink(*[f'{input_path}.{x}' for x in ['bed','bim','fam']])
    elif input_type=='bgen':
        hl.index_bgen(input_path, contig_recoding={"01": "1"}, reference_genome='GRCh37')
        mt = hl.import_bgen(input_path, 
            entry_fields = ['GT', 'GP'], 
            sample_file = '/gpfs1/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample')
    if input_type!='mt' and cache:
        mt = mt.cache()
    return mt


def recalc_info(mt, maf=None, info_field='info', gt_field='GT'):
    r'''Recalculate INFO fields AF, AC, AN and keep sites with minor allele
    frequency > `maf The fields AF and AC are made into integers instead of 
    an array, only showing the AF and AC for the alt allele
    '''
    if info_field in mt.row:
        mt = mt.annotate_rows(
            info = mt[info_field].annotate(**hl.agg.call_stats(mt[gt_field],
                                                               mt.alleles)
                                               ).drop('homozygote_count')
            ) # recalculate AF, AC, AN, drop homozygote_count
    else:
        mt = mt.annotate_rows(**{info_field: hl.agg.call_stats(mt[gt_field], mt.alleles)})
    mt = mt.annotate_rows(
        **{info_field: mt[info_field].annotate(
            **{field:mt[info_field][field][1] for field in ['AF','AC']}
            )}
        ) # keep only the 2nd entries in AF and AC, which correspond to the alt alleles
    if maf is not None:
        mt = mt.filter_rows((mt[info_field].AF>maf) & (mt[info_field].AF<(1-maf)))
    return mt

def get_fam(app_id=11867, wes_200k_only=False):
    assert app_id in {11867,12788}
    fam_path = f'/well/lindgren/UKBIOBANK/nbaya/resources/ukb{app_id}_{"wes_200k_" if wes_200k_only else ""}pedigree.fam'
    fam = hl.import_table(paths=fam_path, key='IID',types={f:'str' for f in ['FID','IID','PAT','MAT','SEX','PHEN']})
    return fam

def filter_to_unrelated(mt, get_related=False, maf=None):
    r'''Filter to samples in duos/trios, as defined by fam file
    '''
    fam = get_fam()
    fam = fam.filter(hl.is_defined(mt.cols()[fam.IID])) # need to filter to samples present in mt first before collecting FIDs
    fam_ct = fam.count()
    col_ct = mt.count_cols()
    print(f'\n{col_ct-fam_ct}/{col_ct} samples have IIDs missing from fam file') # samples with missing IIDs cannot be included in the related sample set
    iids = fam.IID.collect()
    offspring_fids = fam.filter(hl.literal(iids).contains(fam.PAT)
                                |hl.literal(iids).contains(fam.MAT)).FID.collect() # get the FIDs of offspring with at least one parent in the dataset
    fam = fam.filter(hl.literal(offspring_fids).contains(fam.FID)) # subset to samples which share FIDs with offspring
    if get_related:
        mt = mt.filter_cols(hl.is_defined(fam[mt.s])) # related
    else:
        mt = mt.filter_cols(~hl.is_defined(fam[mt.s])) # unrelated
    mt = recalc_info(mt=mt, maf=maf)
    return mt

def incorporate_phenotypes(mt, pheno_path):
    ht = hl.import_table(paths=pheno_path, delimiter =',', key = 'ID', impute = True, types = {'ID': hl.tstr})
    mt = mt.annotate_cols(pheno=ht.index(mt.s))
    print(f'\n{pheno_path} was sucessfully loaded..')
    return mt

def subset_variant(mt, variant):
    v = hl.parse_variant(variant)
    mt = mt.filter_rows( (mt.locus == v.locus) & (mt.alleles == v.alleles) )
    #mt = mt.filter_rows(mt.locus == v.locus)
    #mt = mt.filter_rows(mt.alleles == v.alleles)
    count = mt.count()
    print(f'\n{count[0]} rows were found for variants "{variant}".')

def linear_regression(mt, phenotype):
    result = hl.linear_regression_rows(
        y=[mt.pheno[phenotype]], 
        x=mt.GT.n_alt_alleles(), 
        covariates=[1, mt.pheno.age, mt.pheno.genotyping_array, mt.pheno.ukbb_centre, mt.pheno.sex,
                       mt.pheno.PC1, mt.pheno.PC2, mt.pheno.PC3, mt.pheno.PC4,
                       mt.pheno.PC5, mt.pheno.PC6, mt.pheno.PC7, mt.pheno.PC8,
                       mt.pheno.PC9, mt.pheno.PC10])
    return result

def main(args):

    # load inputs
    input_path = str(args.input_path)
    input_type = str(args.input_type)
    pheno_path = str(args.pheno_path)
    out_prefix = str(args.out_prefix)
    chrom = int(args.chrom)
    pheno = args.pheno 
    print(pheno)

    get_unrelated = args.get_unrelated
    get_wb = args.get_wb

    # init hail and get phenotypes
    hail_init(chrom)
    mt = get_table(input_path, input_type) # 'data/imputed/ukb_imp_chr17_v3.bgen'
    mt = incorporate_phenotypes(mt, pheno_path)

    # subset samples
    if get_unrelated:
        mt = filter_to_unrelated(mt)
    if get_wb:
        mt = mt.filter_cols(mt.pheno.white_british == 1)
    
    # verbose
    counts = mt.count()
    print(f'\nRunning GWAS on {counts[0]} variants and {counts[1]} samples..')

    # perform linear regression
    for phenotype in pheno:
        print(phenotype)
        result = linear_regression(mt, phenotype) # 'Hand_grip_strength_(left)_combined_white_ritish_InvNorm' 
        result.export(f'{out_prefix}.{phenotype}.tsv')
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--pheno_path', default=None, help='Path to phenotype input')
    parser.add_argument('--pheno', nargs='+', help='What phenotypes should be used?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--chrom', default=None, help='chromsome')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated')
    parser.add_argument('--get_wb', action='store_true', help='Get self-reported white british')
    args = parser.parse_args()

    main(args)





