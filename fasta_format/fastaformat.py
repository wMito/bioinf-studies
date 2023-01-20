import argparse
import re



class FASTAFormat():  
    """
    Handles FASTA files. 
    Splits and saves sequences in a dictionary that can be used to create further summaries. 
    """
    
    
    def __init__(self):
        self.source_file = None
        self.dict = None
    
    def read_FASTA(self, filename, del_gaps):
        """
        Creates self.dict which contains a dictionary with all sequences extracted from the FASTA file.
        It also saves the file path for __repr__ and debugging purposes.
        """
        self.source_file = filename
        with open(filename) as dump:
            f = dump.read()
        
        self.dict = {d:s for (d,s) in zip(self.get_descr(f), self.get_seqs(f, del_gaps))}
    
    
    def __repr__(self):
        return f"FASTA_format('{self.source_file}')"
    
    def __str__(self):
        return f"Found {len(self.dict)} sequences.\n" +   \
               f"Sequences IDs: {', '.join(self.dict.keys())}"
    
    def get_descr(self, fasta):
        """
        Returns a list containing all descriptions from the string.
        """
        descr = re.findall('>.*', fasta)
        return [i[1:] for i in descr]    # Deleting > 

    def get_seqs(self, fasta, delgaps):
        """
        Returns a list containing all sequences from the string.
        """
        seqs = re.split(">.*\n", fasta)[1:]
        seqs = [i.replace('\n', '') for i in seqs]
        if delgaps == True:    # Deleting gaps
            seqs = [i.replace('-', '') for i in seqs]
        return seqs 
    
    
    def get_comp_aa(self):
        """
        Creates a dictionary comp_aa containing seq IDs and dictionaries with amino acid compositions.
        """
        self.comp_aa = {}
        all_aa = list("ACDEFGHIKLMNPQRSTVWY")
        
        for sid in self.dict:
            seq_len = len(self.dict.get(sid).replace('-', ''))
            comp_dict = {aa : round(self.dict.get(sid).count(aa)/seq_len, 4) 
                         for aa in all_aa}
            
            self.comp_aa[sid] = comp_dict
      
    def str_comp_aa(self, sid):
        """
        Returns a human readable string summarising amino acid composition for the selected sequence.
        """
        summary = [f"{i}: {self.comp_aa[sid].get(i)}" for i in self.comp_aa[sid]]
        return "Chosen sequence ID: " + sid + "\n" +   \
                "The amino acid composition for this sequence: \n" + '\n'.join(summary)
                
                 
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='input FASTA file')
    parser.add_argument('-g', '--delgaps', help='should gaps be deleted', action='store_true')
    args = parser.parse_args()
    
    fasta = FASTA_format()
    fasta.read_FASTA(args.file, args.delgaps)

    first_seq = list(fasta.dict.keys())[0]
    fasta.get_comp_aa()
    print(fasta.str_comp_aa(first_seq))
    
if __name__ == "__main__":
    main()