# CherSNP

Some useful URL:
https://triticeaetoolbox.org/wheat/tutorials/Tutorial_Variant_Effect.pdf
http://www.sequenceontology.org/browser/current_release/term/SO:0001819

About 0,1 based
https://www.biostars.org/p/84686/
BED is 0-based, VCF is 1-based, GFF is 1-based

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4128822/
Some ideas for training

https://www.nature.com/articles/s41587-019-0164-5
CUG and GUG start codon...

The consensus sequence for an intron (in IUPAC nucleic acid notation) is: G-G-[cut]-G-U-R-A-G-U (donor site) ... intron sequence ... Y-U-R-A-C (branch sequence 20-50 nucleotides upstream of acceptor site) ... Y-rich-N-C-A-G-[cut]-G (acceptor site).


		Some problem...:
			1. In GRCh38, some predicted transcripts only have 1 CDS and nothing else
				which would break the training method
				Probably just skip all XM_***?

https://www.sciencedirect.com/science/article/pii/S0022283697909517?via%3Dihub

http://genesdev.cshlp.org/content/31/17/1717.full

python3 tester.py --fasta data/foo.fasta --dict default/Jeanne.js --init default/Katyusha.pkl --term default/Erika.pkl --entCDS default/Nadeshiko.pkl --outCDS default/Juliet.pkl

# One nucelotiude start codon
NM_001110358.1 AGT [[2522, 2522], [359, 990]]
NM_001110360.1 AGT [[2447, 2447], [359, 990]]
NM_001110361.1 AGT [[2522, 2522], [359, 990]]
NM_001110359.1 AGT [[2447, 2447], [359, 990]]

Introduction in short:

My project is aiming to come out a program which can predict the effects of SNPs within given transcription sequence like Ensemble's VEP could do but neither require users to download large data cliche nor continuously interact with server.

To achieve this goal, unlike VEP and other widely known prediction's method of finding appropriate reference sequences and then using documented annotations, my program would apply SVM technique to determine important sites within the given sequence and use predicted positions of the sites to generate correspond annotations. After that, my program will stimulate translation process to eventually predict the effects on residue level.

So far, I have developed prototype methods to train classifiers to identify candidate start sites, stop site, donor site, and accepter site which can help my program to identify possible CDSs within the input sequence. Then, correspond score of each candidates will be calculated and candidates with score no less than first quantile will be considered as high quality candidates.

By now, I am still tying to find out the best threshold to get high quality candidates.

Once the high quality candidates acquired, my program would check out the possible combinations of start sites, stop sites, and introns, generate possible CDS patterns, and then eliminate which does not make senese biologically, (e.g. multiple stop codon in reading frame, CDS total length not in order of 3...)

With validated CDS pattern, my program will start stimulating translation process and determine correspond residue changes of given SNPs and output prediction in SO terminology which consist with what VEP would output.

Since last update, the basic framework has been set up with agreement from both professor Lessick and Professor Yun:

When training SVM classifiers, either sequences labeled predicated (e.g. have name start with "X") or sequences having CDS region with length less than 3 will be disregarded. 4 classifiers will be trained to out to identify start site, donor site, acceptor site, and stop site accordingly.

After classifiers are trained out and identified potential sites within query sequences, potential transcript patterns will be calculated out through testing every combinations of annotations and the patterns which do not make sense with transcription manner (e.g. Overall coding sequence is not in factor of 3 and in frame stop codon appears before stop site) will be disregarded. With potential transcript patterns, translation process will be stimulated and the effect of SNP will be predicted based on whether or how the translation process would be interrupted.

With suggestion from Professor Yun, the niche of my project will focus more on predicting effects of variants in sequences which do not have generally faithful annotations. Based on this idea, my program will be later tested with Cov19 variants and at this moment I am still studying about how variants could affect the normal function of virus. and maybe also bacteria later depends on time available.
