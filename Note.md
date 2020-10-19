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
