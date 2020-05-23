# CherSNP
Predict Effect of SNPs based on HGVS annotation and DNA sequence (or transcript sequence)
without having internet access nor large cache data

Find CDS within Sequence
SVM train to find splite sites?
Would depend on Sci-kit


Some useful URL:
https://triticeaetoolbox.org/wheat/tutorials/Tutorial_Variant_Effect.pdf
http://www.sequenceontology.org/browser/current_release/term/SO:0001819

About 0,1 based
https://www.biostars.org/p/84686/
BED is 0-based, VCF is 1-based, GFF is 1-based

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4128822/
Some ideas for training

There would be 4 seperate classifier under 2 major type:
	Type-in: (Aoko)
		1 Start of cds, should have initial codon down stream
		2 End of intron, go back to exon

	Type-out: (Alice)
		1 End of exon, going in to intron
		2 End ot all cds, should have terminal signal upstream

	Aoko takes 20 bp upstream and 6 bp downstream:
		1 	should have 5'UTR containing protein binding sites (GC rich expected?)
			Start codon should some how looks like AUG/ATG
			What codon comes after start codon might also be a good stuff?
		2	End of intron has CT(Y)rich region for splice
			Intron usually ends with AG also should have something in exon more ofenly

	Alice takes 6 bp upstream and 4 bp downstream:
		1	Intron tends to start with GU and also something in exon for splicing purpose
		2	Terminal codons are also pretty constrained, +4 bp downstream context and codon 		before terminal codon would be helpful?

	What dimention classifier would have?
		Aoko1: GC ratio in upstream, start codon, second codon
		Aoko2: CT ration in upstream, end of intron, 4bp downstream context
		Alice1: n bp upsteam, start of intron, 2bp downstream
		Alice2: -2 codon, terminal codon, 4bp downstream context

	How to train?
		Again, make a dictionary of motif frequency for both splice site and non splice sites.
		Compare frequecy difference within 2 dictionaries, try to find the "smallest big probability of motif which can determine 2 types of sequence", and set these values for SVM kernel (Linear?) or use RBF and polynomial on scikit?
