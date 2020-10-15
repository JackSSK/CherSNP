CherSNP Project Summary (in progress)

Objective:
  Aiming to come out with a package being able to predict the effect of Single-nucleotide polymorphisms(SNPs) by applying machine learning methods.

Introduction:
 	When I was doing an internship in a research institution focus on building Autism related gene database few months ago, I participated in developing programs for standardize genomic variant annotations and put them into new variant database. One part I was asked to do is predicting the effects of variants and residue changes (change in protein level) for which only have allele change annotations (change in DNA/RNA level).  
	Firstly, I was suggested to use a R package, VariantAnnotation, to implement; however, due to the very restricted prediction output types (only synonymous or non-synonymous) which do not meet the nomenclature standard of database, we eventually chose to use Ensembl’s Variant Effect Predictor (VEP) instead.

  Even though VEP certainly can give very decent and accurate predictions, in my opinion, there are still some inconvenient aspects of using it in a scenario of predicting effect of large amount of variants in few genes:
  	1. If using VEP in online mode, the running time for predictions could be unstable and, depends on connecting status, the running time could be unacceptably long.
  	2. If using VEP in offline mode, the running time is indeed much stable and clearly shorter, but this method requires users to download large datasets with size depends on how many reference genomes needed. (In case of human  GRCh37 and GRCh38 reference genomes, I needed to downloaded 25+ Gb datasets.)

  Another aspect of VEP could be further improved is that the accuracy and reliability of the predictions could vary significantly depends on sequence being picked as reference. For example, if the prediction is made based on a sequence labeled Matched Annotation from NCBI and EMBL-EBI (MANE), it's accuracy and reliability would be much more valued than prediction made based on predicted sequences. In the case of all predictions of a variant are made based on the predicted sequences, the quality of predictions could be questionable. After I checked for SnpEff and ANNOVAR which I think could be considered as alternative choices of VEP, I found similar inconveniences also exist in them. Therefore, I am wondering whether I could develop a human variant effect predictor through using machine learning methods to identify locations of features such as Coding Regions(CDS), Introns, start sites, and stop sites. With such method, neither connection with server nor downloading large datasets would be required, and it would be able to output predictions having relatively consistent qualities. Furthermore, pilot studies could be done more easily, and the relatively low quality predictions made by VEP could be double checked.

  By now, this project only apply Support Vector Machine (SVM) method to train for classifiers. Later, in order to achieve higher accuracy, Neural Network (NN) methods ,which is also implemented in the package, sklearn I used, and Hidden Markov Model methods may also be applied.

Method:
	As a result of not asking users to download large datasets, the predictor would need users to provide reference sequences along with variant annotations at least. Thus, the minimum input needed should be a FASTA file which have reference sequences and a GFF/BED file which has annotations in HGVS format [1]. Ideally, user also includes CDS features in the GFF/BED, but this is not mandatory and, in this case, trained classifiers will be used to identify CDSs, introns, Start Sites, Stop Sites within the input sequences.

  To train out such classifiers, users may provide sample sequences with associated annotations for CDS locations in FASTA and GFF (maybe also bed later) format accordingly. With examples, four types of classifiers would be trained in order to locate the start sites, stop sites, and boundaries between CDSs and introns:

    1. Type init: Classifier to identify where the sequence going from UTG into CDS, the location of start site. The concept of 'Kozak consensus sequence' will be applied which could be helpful for not missing start sites not having start codon being ATG/AUG. 4bp upstream, start codon, and 3bp downstream of documented start sites (area around first codon in first CDS) will be used as positive samples and 2 sequences generated through forward or backward shifted by 3bps as negative samples.

    2. Type outCDS: Classifier to identify where the sequence going from CDS into intron. Considering consensus structure of introns, the 2bp upstream of cutting position (which should be like GG [3]), first 2bp in intron (should be like GT [3]), and next 4bp in intron (likely being consensus with RAGU [3]) will be documented as positive sample of donor site, and 2 negative samples will also be gained through 3bp shift to either side.

    3. Type entCDS: Classifier to identify where the sequence going from intron into CDS.
    Considering consensus structure of introns, the Y-ratio of 15bp to 2bp upstream of cutting position (expected to be high[3]), 2bp upstream (should be like AG[3]), and 1bp in adjacent CDS will be documented as positive sample of acceptor site, and 2 negative samples will be gained through 3 bp shifts as well.

    4. Type term: Classifier to identify where the sequence going from CDS into UTG, the location of stop site. Now the classifier training still process like Type init, which take 3bp upstream (last Amino Acid), stop codon, and 4bps downstream at area around last codon in last CDS. Positive and negative samples are also gained from similar manner. However, this design is mainly due to considerations of users studying species which may have stop codon also coding for amino acid like ciliates.[4] For human studies, which would be the main purpose as I design it, this may not be necessary and even actually lower accuracy than simply detecting for stop codon. Thus, modifications may be made later.

  After positive samples being stored in a hash table (dictionary in python), each type of subsample's associated number of times being observed from training data would be transposed into log scale through:

          log(# of time / total observation)

  The table, dict, gained after transposition will be used to transform read-in sequences into numbers not only in training process but also in predicting process, and it will be stored into an output file after training.

  By using dict, positive samples and negative samples will be fit into 4 SVM classifiers through using rbf kernel by now. (Later may also try out other methods like NN and HMM) and the resulted classifiers will be stored into 4 separate files.

  With classifiers, the prediction part will start with identifying possible locations of start sites, cdss, introns, and stop sites accordingly, and list of possible alternative splicing methods will be gained through combining features based on their positions. (e.g. start sites must go first) and not contradicting with typical gene/transcript structure (e.g. start site cannot overlap with donner site).

  Then, the SNP annotations will be decoded and applied to list of candidates. The associated effects will be gain through stimulating translation process with each candidate's features positions. The nomenclatures of predicted effects will be consisting with MISO[2], which is consisting with VEP's nomenclatures.



Reference:
[1] https://varnomen.hgvs.org/
[2] http://www.sequenceontology.org/browser/current_release/term/SO:0001819
For prediction type: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
[3] "Molecular Biology of the Cell". 2012 Journal Citation Reports. Web of Science (Science ed.). Thomson Reuters. 2013.
[4]https://www.cell.com/cell/pdf/S0092-8674(16)30788-7.pdf
