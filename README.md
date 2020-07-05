Jack Yu's thesis project proposal



Introduction:
 	When I was doing an internship in a small research institution focus on building Autism related gene database few months ago, I participated in developing programs for standardize genomic variant annotations and put them into new variant database. One part I was asked to do is predicting the effects of variants and residue changes (change in protein level) for which only have allele change annotations (change in DNA/RNA level).  
	Firstly, I was suggested to use a R package, VariantAnnotation, to implement; however, due to the very restricted prediction output types (only synonymous or non-synonymous) which do not meet the nomenclature standard of database, we eventually chose to use Ensembl’s Variant Effect Predictor (VEP) instead. Even though VEP certainly can give very decent prediction, I think there are still some inconvenient aspects of using it in a scenario of predicting effect of large amount of variants in few genes:
	1. If we use VEP in online mode, the running time for prediction was kind of unacceptable.
	2. If we use VEP in offline mode, the running time is indeed much shorter, but we needed to downloaded 25+ Gb datasets for GRCh37 and GRCh38 reference genomes.
	After I checked for SnpEff and ANNOVAR which I think could be considered as alternative choices of VEP, I found similar inconveniences also exist in them. Therefore, I am wondering whether I could develop a human variant effect predictor in a different approach which would neither require need users to connect with online database nor download large datasets.



Method:
	As a result of not asking users to download large datasets, the predictor would need users to provide reference sequences along with variant annotations at least. Thus, the minimum input needed should be a FASTA file which have reference sequences and a GFF/BED file which has annotations in HGVS format [1]. Ideally, user also includes CDS features in the GFF/BED, but this might be hard for users to do and I would not make it mandatory.

	Here is my imagination on prediction making process so far:
	1. Process GFF/BED file first and then save variant information in arrays based on reference sequences. If CDS annotations are available, also save them.
	2. Process sequences in FASTA one by one. If CDS annotations are available, use so. Otherwise, use trained classifiers to label out the CDSs. Then, label out 5’UTR, 3’UTR, 	introns, splice acceptors, and splice donors according to CDS coordinates. Also, stimulate translation and save the result.
	3. For each variant of the reference sequence, check the position of variant and adjacent base pairs and then make prediction based on translation stimulating manner. The predicted effect should be in nomenclature same with MISO. [2]

	And here is my imagination on classifier training:
	1. In total, there should be 4 types of classifiers which can determine 4 types of changes:
		5’ UTR → CDS (Type UC)
		CDS → Intron (Type CI)
		Intron → CDS (Type IC)
		CDS → 3’UTR (Type CU)

	2. Even though so far I am planning to Support Vector Machine (SVM) to train for all 4 types of classifiers, different classifiers would be train with different features.
		Type UC, I would use GC ratio in a selective region before the cutting edge and 	possibility of having the first codon after cutting edge (because the TF binding sites before 	would have a high GC ratio and start codon should be AUG or something close)
		Type CI, I would use the possibility of having the 2 bps before the cut and 2bps after the 	cut and 4 bps further. (The consensus sequence would be GG/GURAGU, but probably the GU 	feature immediately after cut is more stronger?)
		Type IC, I would use CT ratio in a selective ratio before the cutting edge and possibility 	of having last 2 bp before cutting edge. (intron usually end with AG and have CT rich region 	slightly before)
		Type CU, I would just use the possibility of having the last codon as terminal codon.

By now, I imagine the predictor will eventually be a python package and may relay on Scikit-learn for SVM features.


Reference:
[1] https://varnomen.hgvs.org/
[2] http://www.sequenceontology.org/browser/current_release/term/SO:0001819
For prediction type: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
