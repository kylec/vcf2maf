import requests
import argparse

# annotate maf with oncotater webapi. annotation columns are appeneded to maf file
# annotation columns 
# gene   transcript    variant_class      dbsnp       aa_change  pph2_class
# OR4F5	 NM_001005484  Missense_Mutation  rs75062661  p.T141A	 neutral
#
# Kyle Chang

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input')
    parser.add_argument('-o', dest='output')
    args = parser.parse_args()
    
    url = 'http://www.broadinstitute.org/oncotator/mutation/'
    
    # open files
    try:
        mafFile = open(args.input, 'r')
    except IOError:
        print 'Error: maf input does not exist.'    
    annoFile = open(args.output, 'w')

    # headers
    annoHeader = 'gene\ttranscript\tvariant_class\tdbsnp\taa_change\tpph2_class'
    mafHeader = mafFile.readline().rstrip('\n')
    annoFile.write(mafHeader + '\t' + annoHeader + '\n')
    
    # read maf file for chr, start, end, ref, var
    try:
        i = 0
        for row in mafFile:
            i += 1
            fields = row.rstrip('\n').split('\t')
            query = url + fields[4] + '_' + fields[5] + '_' + fields[6] + '_' + fields[10] + '_' + fields[12]
    
            # query oncotator
            query = requests.get(query)
            result = query.json()
    
            # get best transcript and its annotation
            txId=result['best_canonical_transcript']
            anno=result['transcripts'][txId]
    
            # get polyphen results
            if result.has_key('pph2'):
                pph2=result['pph2']['pph2_class']
            else:
                pph2 = ''
        
            if result.has_key('dbSNP_RS'):
                dbsnp = result['dbSNP_RS']
            else:
                dbsnp = ''
        
            if anno.has_key('refseq_mRNA_id'):
                transcript = anno['refseq_mRNA_id']
            else:
                transcript = ''
        
            # write annotation to output    
            fields.extend([anno['gene'], transcript, anno['variant_classification'], dbsnp, anno['protein_change'], pph2]) 
            annoFile.write('\t'.join(fields))
            annoFile.write('\n')
    except Exception:
        print 'Error in maf record %d' % i
    
    mafFile.close()
    annoFile.close()

if __name__ == '__main__':
    main()