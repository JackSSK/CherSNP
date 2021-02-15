# CherSNP ver.0.2
I broke everything

  Don't try to run all the code now!
  
If you insist:
  python3 train.py --fasta data/Chr14.fasta --gff data/Chr14.gff --mode hasRNA
  python3 tester.py --fasta data/TTC5.fasta
  And it works much better in the case of Covid instead of human genome

Compared with ver.0.1, we need 2 major changes
1st: More generalized frame -Anze
  Have classes and functions only defined like:
    def functionOne (input):
      return output
  As the generalized frame defined, we will populate the classes and functions with codes in ver.0.1
  
2nd: New classifiers
  Someway somehow utilize gapped-kmer method with part of the seq being constrained(if necessary)
  1st question: Do an estimation on required computational power if we use the method:

    一： 
          Given key (if there is one)
          vectors = []
          for feature f:
            for all observed motif m in look_up_dict[f]:
              a = a coefficient associated with how many times m has been observed in training set
              v = # of shared gapped k-mer(m, input) * a 
              vectors.append(v)
    二：
          Given key (if there is one)
          vectors = []
          for feature f:
            m = most commonly observed motif in look_up_dict[f]
            a = a coefficient associated with how many times m has been observed in training set
            v = # of shared gapped k-mer(m, input) * a 
            vectors.append(v)
  
  Tune the parameters little bit and try to find how performances could be changed: -Kano?
    What if we combine both method through only use top 10% of the commonly seen motifs
    What if we have larger k values like k >= 6
    What if we have each motif m being longer like len(m)>=10
