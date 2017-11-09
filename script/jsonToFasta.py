import json
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict

def jsonToFasta(input_file):
    data = None
    with open(input_file) as handler:
        data = json.load(handler)
        handler.close()
    new_records = []
    with open(input_file) as handler:
        data = json.load(handler)
        handler.close()
    for serotype, alleles in data.items():
        for allele in alleles:
            new_record = SeqRecord(Seq(allele['seq'],IUPAC.IUPACUnambiguousDNA),
                                   id=serotype+'-'+str(allele['num']),
                                   description="")
            new_records.append(new_record)
    SeqIO.write(new_records, "Data/serotype_dict.fasta", 'fasta')

def mimicDictionary(input_file):
    data = None
    with open(input_file) as handler:
        data = json.load(handler)
        handler.close()
    new_dict = defaultdict(dict)
    for serotype, alleles in data.items():
        for allele in alleles:
            serotype_class = serotype[0]
            new_entry_name = serotype+'-'+str(allele['num'])
            new_entry = {"allele": serotype, "gene": allele['gene']}
            new_dict[serotype_class][new_entry_name] = new_entry
    with open('Data/allele_serotype.json', 'w') as handler:
        json.dump(new_dict, handler, indent=4, separators=(',', ': '))
        handler.close()


def main():
    jsonToFasta('Data/serotype_dict.json')
    mimicDictionary('Data/serotype_dict.json')

if __name__ == '__main__':
    main()