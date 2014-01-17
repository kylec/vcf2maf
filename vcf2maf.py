import vcf
import argparse

# read snp entries in vcf, output a maf
# depend on PyVCF. pip install PyVCF
# Kyle Chang

# define maf headers
mafFields = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',                             
    'Chromosome', 'Start_position', 'End_position', 'Strand',     
    'Variant_Classification', 'Variant_Type', 'Reference_Allele',
    'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 
    'Tumor_Sample_Barcode', 'Normal_Norm_SAmple_Barcode', 
    'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
    'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 
    'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 
    'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID']
#TODO add annotation info

# define variant status for MAF    
variantStatus = {0 : 'Wildtype', 1: 'Germline', 2 : 'Somatic', 3 : 'LOH', 4 : 'post-transcriptional modification', 5 : 'Unknown' }

# convert vct genotype(0/0,0/1,1/1) to maf alleles (A T)
# homo var when zero 0, i.e. 1/1
# het when one 0, i.e. 0/1
# homo ref two 0, i.e. 0/0
def toMafAlleles(vcfGenotype, refAllele, varAllele):
    count = vcfGenotype.split('/').count('0')
    
    if count == 0:
        allele1 = allele2 = varAllele
    elif count == 1:
        allele1 = refAllele
        allele2 = varAllele
    else:
        allele1 = allele2 = refAllele
    
    return (allele1, allele2)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='vcf')
    parser.add_argument('-o', dest='maf')
    args = parser.parse_args()
    
    # init maf fields
    tumorBarcode = normalBarcode = tumorUuid = normalUuid = platform =  '.'
    tumorAllele1 = tumorAllele2 = normalAllele1 = normalAllele2 = '.'
    
    vcfFile = open(args.vcf, 'r')  
    vcfReader = vcf.Reader(vcfFile)
    mafFile = open(args.maf, 'w')
    
    # vcf metadata = ordered dict
    # metadata['center'] = list
    center = vcfReader.metadata['center'][0]
    #ncbiBuild = vcfReader.metadata['reference']
    
    # returns a list of samples metadata, each sample metadata is an ordered dict 
    samples = vcfReader.metadata['SAMPLE']
    for sample in samples:
        
        if sample['Platform']:
            platform = sample['Platform'] 
        
        if 'normal' in sample['Description'].lower():
            normalBarcode = sample['SampleTCGABarcode']
            normalId = sample['ID']
        elif 'rna' in sample['Description'].lower():
            # skip RNA tumor barcode, seen in UCSC VCF
            print "Skip RNA tumor..."
        elif 'tumor' in sample['Description'].lower():
            tumorBarcode = sample['SampleTCGABarcode']
            tumorId = sample['ID']
    
    # read vcf records
    for vcfRow in vcfReader:
        chrom = vcfRow.CHROM
        start = str(vcfRow.POS)
        end = str(vcfRow.POS)
        
        # variant type (SNP, DEL, DNP etc)
        # TODO allow all variant types
        varType = vcfRow.INFO['VT']
        if varType != 'SNP':
            continue
        
        # reference, variant, tumor, and normal alleles
        ref = vcfRow.REF
        var = str(vcfRow.ALT[0])
        normalAllele1, normalAllele2 = toMafAlleles(vcfRow.genotype(normalId)['GT'], ref, var)
        tumorAllele1, tumorAllele2 = toMafAlleles(vcfRow.genotype(tumorId)['GT'], ref, var)
       
        # Mutation Status (Somatic/Gemrline/LOH/Wildtype)
        mutationStatus = variantStatus[vcfRow.genotype(tumorId)['SS']]
        
        # write MAF
        mafRow = []
        mafRow.append('')               #hugo_symbol
        mafRow.append('')               #entrez_gene_id
        mafRow.append(center)
        mafRow.append('')               #mcbi_build
        mafRow.append(chrom)
        mafRow.append(start)
        mafRow.append(end)
        mafRow.append('+')
        mafRow.append('')               #variant_classification
        mafRow.append(varType)
        mafRow.append(ref)
        mafRow.append(tumorAllele1)
        mafRow.append(tumorAllele2)
        mafRow.append('')               #dbsnp
        mafRow.append('')               #dbsnp_valstat
        mafRow.append(tumorBarcode)
        mafRow.append(normalBarcode)
        mafRow.append(normalAllele1)
        mafRow.append(normalAllele2)
        mafRow.append('')               #validaiton alleles
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')               #validation_status
        mafRow.append(mutationStatus)
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')
        mafRow.append('')
        mafRow.append(platform)
        mafRow.append('')
        mafRow.append('')
        
        for field in mafRow:
            #print field, type(field)
            mafFile.write(field)
            mafFile.write('\t')
        mafFile.write('\n')
    
    vcfFile.close()
    mafFile.close()
if __name__ == '__main__':
    main()
    