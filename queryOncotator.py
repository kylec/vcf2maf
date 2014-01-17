import requests

# test query to onocator

r = requests.get('http://www.broadinstitute.org/oncotator/mutation/chr7_140453136_140453136_A_T')
anno = r.json()
txId=anno['best_canonical_transcript']
tx=anno['transcripts'][txId]
print tx['protein_change']