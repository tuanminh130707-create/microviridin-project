# Identifying Functionally Important Residues in Graspetide ATP-Grasp Macrocyclases (MdnB and MdnC): A Literature Review

*This review pulls together two bodies of work: (1) what is known about the biology and structure of the microviridin/graspetide macrocyclases MdnB and MdnC, and (2) the established computational methods for pinpointing which amino-acid positions in an enzyme matter for its function. The aim is to justify, from the literature, a defensible workflow for finding functionally important residues in these enzymes. Technical terms are defined on first use and collected in a **Glossary** at the end.*

---

## 1. Introduction

Ribosomally synthesized and post-translationally modified peptides — **RiPPs** — are a large class of natural products built by the cell's normal protein-making machinery and then chemically remodeled by dedicated tailoring enzymes. Among them are the **graspetides** (historically called microviridins), small peptides locked into rigid cages by internal cross-links, many of which are potent inhibitors of serine proteases and therefore of interest as drug leads.

The cross-links that define a graspetide are installed by **macrocyclases** — enzymes belonging to the **ATP-grasp ligase** family. In the founding microviridin system there are two of them: **MdnC**, which forms two ester cross-links, and **MdnB**, which forms one amide cross-link. Understanding *how* these enzymes recognize their substrate and catalyze these reactions comes down to a smaller question: which of the ~320 amino-acid positions (**residues**) in each enzyme actually matter, and for what? Those positions are the targets a researcher would mutate to test mechanism or to re-engineer the enzyme.

This review surveys (Section 2) the structural biology of MdnB and MdnC and their close relatives, (Section 3) the genome-mining work that has revealed how many related sequences exist to learn from, and (Section 4) the computational methodology for turning those sequences into a ranked list of important residues, before synthesizing these into a workflow (Section 5).

## 2. The Biosynthetic System and Its Enzymes

### 2.1 Precursor logic: leader and core

Graspetides begin life as a longer **precursor peptide** (in microviridins, MdnA). This precursor has two functionally distinct parts. The **core peptide** is the segment that ends up in the final product and carries the residues that become cross-linked. The **leader peptide** is an N-terminal segment that is *not* part of the final product but is essential for recognition: the enzymes bind the leader, not the core. In microviridins the leader carries a strictly conserved motif, **PFFARFL**, that folds into a short helix and docks onto the enzyme (Li et al., 2016). This "leader-directed" logic — recognition and catalysis handled by different parts of the substrate — is a recurring theme across RiPP biosynthesis.

### 2.2 The two macrocyclases and their chemistry

MdnC and MdnB are homologous ATP-grasp ligases that act in sequence. **MdnC installs the two ester cross-links** (technically **macrolactones**, or ω-esters, formed between a Ser/Thr hydroxyl and an Asp/Glu carboxyl side chain); **MdnB then installs the single amide cross-link** (a **macrolactam**, or ω-amide, between a Lys amine and an Asp/Glu carboxyl) (Li et al., 2016; Ramesh et al., 2021). (In a second, equally common naming scheme these enzymes appear as **MvdD** = the lactone-former and **MvdC** = the lactam-former; the literature uses both sets of names for the same two functions.)

Mechanistically, both enzymes use **ATP** and **Mg²⁺**. The ATP-grasp fold "grasps" a molecule of ATP and uses it to phosphorylate the acceptor carboxylate, generating a reactive **acylphosphate intermediate** that the incoming nucleophile (a hydroxyl or amine) then attacks to form the cross-link, releasing ADP and phosphate (Ramesh et al., 2021). This same chemistry is used across the wider ATP-grasp family by primary-metabolism enzymes such as D-Ala–D-Ala ligase and glutathione synthetase — a point that matters later, because it means the ATP-handling residues are shared across a vast superfamily and are therefore *not* what makes these enzymes graspetide-specific.

### 2.3 Structural insights: the MdnB and MdnC crystal structures

The defining structural study of this system solved crystal structures of both enzymes: MdnB alone and MdnC bound to the MdnA leader peptide (Li et al., 2016). These are deposited in the **Protein Data Bank (PDB)** as **5IG8 (MdnB)** and **5IG9 (MdnC with MdnA)** and are the templates onto which any residue-importance analysis of this system should be mapped. Several findings from that work anchor the rest of this review:

- **Leader recognition is electrostatic and helix-mediated.** The MdnA leader helix contacts a cluster of residues on a helix of MdnC (referred to as **α7**): Glu191, Asp192, and Asn195, which engage Arg17 of the leader. Reversing the charge of two of these (the E191K/D192K double mutant) abolished both binding and catalysis, confirming that this small surface patch is the critical recognition element.
- **A large motion gates the active site.** A structural element called the **β9β10 hairpin** shifts by ~25 Å between the two enzymes: MdnC (bound to leader) adopts an "open" state that admits the substrate, whereas MdnB was captured in a "closed" state. This conformational change, rather than a static pocket, is central to how the enzymes work.
- **The catalytic/ATP residues were assigned by structural comparison.** Because no ATP or nucleotide could be captured in either crystal, the ATP- and metal-coordinating residues of MdnC (Lys125, Lys166, Gln207, Glu215, Asp281, Glu294, Asn296) were identified by aligning the structure with other ATP-grasp enzymes. These should be treated as high-confidence but model-derived assignments.

### 2.4 Comparative structures across the graspetide family

Since 2016, structures of several related graspetide macrocyclases have refined the picture of what is conserved versus what diverges:

- **PsnB** (from plesiocin biosynthesis) was captured with its core substrate and nucleotide, showing that a conserved active-site arginine (**Arg213**) recognizes the ring-forming acidic residue and delivers it toward ATP (Song et al., 2021).
- **CdnC** (chryseoviridin) was solved as a four-part complex with ADP, leader, and substrate peptide, reinforcing the leader-recognition mode and the substrate-induced domain motion (Zhao et al., 2021).
- **PruB** (prunipeptin), a single cyclase that makes both ester and amide links, established through structure and mutagenesis that a conserved **DxR motif** (an Asp–any–Arg pattern; the arginine is PruB Arg233) is critical for catalysis (Rubin et al., 2024). This arginine is the family-wide analog of PsnB's Arg213.
- **Cooperative multi-enzyme systems** further illustrate the range: a "distributive" cyclase processes several core repeats on one long precursor (Zhang et al., 2018), and a recently described deep-sea system uses four ATP-grasp ligases acting cooperatively to build a pentacyclic product, with some ligases active only in the presence of their partners (Li et al., 2026).

Taken together, these structures predict what a conservation analysis of MdnB or MdnC *should* recover: a deeply conserved ATP-grasp catalytic core (shared with the whole superfamily), a family-specific conserved patch on the leader-binding helix (the Glu191/Asp192-type cluster), and a conserved arginine in the acceptor-recognition cleft (the DxR/Arg213-type residue). Importantly, the residues that *distinguish* the lactone-forming (MdnC) clade from the lactam-forming (MdnB) clade are expected to concentrate in the mobile β9β10/α7 elements and the substrate-selecting pocket — the divergent evolution of these two enzyme types being a documented feature of the family (Lee et al., 2020).

## 3. Where the Sequences Come From: Genome Mining and Family Expansion

A conservation analysis is only as good as the set of related sequences (**homologs**) feeding it, so it is worth understanding how the field assembled the large sequence collections now available.

Standard similarity search (**BLAST**) struggles to find RiPP precursors, because the small precursor peptides are short, hypervariable, and frequently not annotated as genes at all. The **RODEO** tool was introduced to solve this by scoring the conserved *neighborhood* around a precursor — combining profile **HMM** searches (statistical models of a protein family) with **machine learning** — first for lasso peptides (Tietz et al., 2017) and then, in a re-trained form, for lanthipeptides across more than 100,000 genomes (Walker et al., 2020). Those studies quantified the "hidden" problem directly, finding that large fractions of true precursors were entirely unannotated in the databases.

For graspetides specifically, two studies matter most. Lee et al. (2020) used homology searches and **sequence similarity networks** (SSNs — graphs that cluster proteins by similarity) to map the ω-ester peptide family into 12 groups, and — crucially for this project — documented that the macrocyclases have undergone **divergent evolution**, sometimes splitting into two specialized enzymes (as in microviridins) and sometimes fusing into one. Ramesh et al. (2021) then applied a dedicated graspetide RODEO module to the ATP-grasp superfamily and roughly doubled the known diversity, identifying 3,923 high-confidence gene clusters and expanding the family from 12 to 24 groups. More recent large-scale efforts continue to enlarge this space, including surveys of human-microbiome genomes that pair computational discovery with chemical synthesis of the predicted peptides (Zhang et al., 2025), and protocol chapters that formalize how to combine RODEO with genomic-enzymology workflows (Hudson & Ramesh, 2025).

The practical upshot for a residue-importance study is twofold: there are now **more than enough** MdnB- and MdnC-like sequences to support a conservation analysis, but they are drawn from a broad and unevenly sampled space, so they must be curated carefully rather than used raw.

## 4. Methods for Identifying Functionally Important Residues

### 4.1 The core principle

The central idea is simple: positions that evolution has kept unchanged across many related proteins are usually the ones that matter for folding, catalysis, or binding, because changes there were selected against. Measuring per-position **conservation** across a good set of homologs therefore provides a proxy for functional importance. The rest of the methodology exists to make that measurement trustworthy.

### 4.2 Collecting homologs

Homologs are gathered by similarity search — **BLAST**, its iterative cousin **PSI-BLAST**, or profile-based **HMMER/jackhmmer**, which are more sensitive to distant relatives. Modern practice (and the ConSurf server's default) searches a pre-clustered database such as **UniRef90**, which already collapses sequences that are >90% identical (Ashkenazy et al., 2016). A recurring lesson across this literature is that **diversity matters more than raw number**: a set of near-identical sequences, however large, carries little information about which positions are truly constrained.

### 4.3 Removing redundancy

Sequence databases are badly skewed by over-studied organisms — thousands of near-identical genomes of the same species — which, if left in, make a handful of clades dominate the analysis and artificially inflate apparent conservation (a "twin bias"). The standard fix is **CD-HIT**, which clusters sequences at a chosen identity threshold and keeps one representative per cluster (Fu et al., 2012). For conservation work, clustering at **90% identity (CD-HIT 0.9)** is the recommended setting; ConSurf's own guidance places the useful window between 80% and 95%, warning that clustering more aggressively than that can strip out the diversity the analysis depends on (Current Protocols, 2021).

### 4.4 The identity window

Two bounds define a usable homolog set. At the **upper end**, sequences more than ~90–95% identical to the query are redundant and are removed (Section 4.3). At the **lower end** lies the **"twilight zone"** of sequence alignment, at roughly **25–35% identity**, below which two sequences may still be genuine relatives but are so diverged that alignment programs can no longer place their residues in correct correspondence — introducing errors that corrupt any downstream conservation score (Rost, 1999). A practical working window is therefore about **35% to 90–95% identity** to the query.

### 4.5 How many homologs

Guidance for **ConSurf**, the most widely used conservation-mapping server, is explicit: at least 5 homologs are required, but meaningful rate estimates need many more, and the recommended range is **50 to 300** sequences. Fewer than ~50 gives statistically noisy estimates; more than ~300 slows the computation with diminishing returns unless the extra sequences add real diversity (Current Protocols, 2021; Ben Chorin et al., 2020).

### 4.6 Aligning the sequences

The homologs are combined into a **multiple sequence alignment (MSA)** — the grid in which equivalent positions across all sequences are placed in the same column. Alignment quality directly determines the reliability of the conservation score, so a high-accuracy aligner is warranted; **MAFFT** is the ConSurf default and a common choice, with MUSCLE and Clustal Omega as alternatives (Ashkenazy et al., 2016).

### 4.7 Scoring conservation

ConSurf computes a per-position **evolutionary rate** using the **Rate4Site** algorithm under an **empirical Bayesian** framework — a statistical method shown to outperform simpler maximum-likelihood estimates, especially for modestly sized datasets (Mayrose et al., 2004). The output is a **conservation grade from 1 to 9**, where **1 = most variable, 5 = intermediate, and 9 = most conserved** (positions with too little data are flagged separately) (Ashkenazy et al., 2016, 2010). These grades can be colored directly onto the 3D structure.

### 4.8 Mapping to structure and choosing targets

The final step maps grades onto the structure (here, 5IG8 for MdnB and 5IG9 for MdnC) and reads off which conserved positions sit where. The key interpretive move for this project is to **distinguish conserved residues in the universal ATP-grasp pocket from conserved residues elsewhere**. The former are conserved across the entire ATP-grasp superfamily and reveal nothing specific to graspetide function; the latter — conserved residues in the **leader-recognition helix and the substrate/acceptor-binding clefts** — are the family-discriminating, mechanistically interesting positions, and thus the most informative mutagenesis targets (consistent with the functional residues mapped by Li et al., 2016; Song et al., 2021; Rubin et al., 2024).

## 5. Toward a Workflow for MdnB and MdnC

Synthesizing the above yields a defensible, literature-grounded pipeline:

1. **Seed** each analysis with the relevant structure's sequence (5IG9 for MdnC, 5IG8 for MdnB) and **collect homologs** by HMMER/PSI-BLAST search of a redundancy-reduced database (UniRef90), using a permissive E-value to reach distant relatives.
2. **Filter** to the ~35–90% identity window (Rost, 1999).
3. **Remove redundancy** with CD-HIT at 0.9 (Fu et al., 2012; Current Protocols, 2021).
4. **Target 50–300 diverse sequences**; relax thresholds if MdnB or MdnC homologs are scarce (Ben Chorin et al., 2020).
5. **Align** with MAFFT.
6. **Score** with ConSurf/Rate4Site and map grades 1–9 onto the structure (Ashkenazy et al., 2016; Mayrose et al., 2004).
7. **Prioritize** conserved (grade 8–9) residues that lie *outside* the ATP pocket — the α7 leader cluster and the DxR/acceptor cleft — as mutagenesis candidates, first confirming the map recovers the already-known functional residues as a sanity check.

Two points specific to this system:

**Analyze MdnB and MdnC separately.** Because the two enzymes form genuinely divergent evolutionary clades (Lee et al., 2020), pooling their homologs into one alignment would average away the very residues that distinguish a lactone-former from a lactam-former. Two independent homolog sets and two ConSurf runs are therefore methodologically necessary, not merely tidier.

**Handling a fused cyclase and in-house sequences (brief).** Where a cyclase is a natural fusion of two ATP-grasp domains, the defensible approach is to split it into its component domains and analyze each against its own homolog set, matching how the field treats these homologous-but-distinct enzymes. Private, in-house sequences can be added to a public homolog set to enrich diversity, but should be treated as ordinary members — a tight cluster of closely related in-house sequences would bias the analysis exactly as an over-sequenced organism does, and CD-HIT clustering (or collapsing them to one representative) mitigates this.

## 6. Conclusion and Outlook

The microviridin/graspetide macrocyclases are unusually well set up for a computational residue-importance study: their structures are solved (5IG8, 5IG9), several functional residues are already experimentally validated, and the sequence family is large and well-characterized thanks to a decade of genome mining. At the same time, systematic conservation-and-structure analyses aimed specifically at prioritizing mutagenesis targets in *individual* graspetide macrocyclases remain thin in the literature — most computational effort in this family has gone into discovery and family classification rather than per-enzyme residue mapping. That gap is the opportunity: a careful, validated conservation map of MdnB and MdnC, benchmarked against the known functional residues and extended to untested positions, would be both a genuine contribution and a reusable template for the broader family.

---

## Glossary

- **Acylphosphate intermediate** — a reactive, phosphate-activated form of a carboxyl group, created using ATP, that the ATP-grasp enzymes form en route to making a cross-link.
- **α7 helix / β9β10 hairpin** — specific structural elements of MdnC; α7 carries the leader-recognition residues, and the β9β10 hairpin is the mobile flap that opens and closes over the active site.
- **ATP-grasp ligase** — an enzyme family whose fold "grasps" ATP and uses it to join two groups by forming a new bond; the graspetide macrocyclases belong to it.
- **BGC (biosynthetic gene cluster)** — a group of neighboring genes that together encode the machinery to make one natural product.
- **BLAST / PSI-BLAST** — tools that search a database for sequences similar to a query; PSI-BLAST iterates to find more distant relatives.
- **CD-HIT** — a program that clusters sequences by similarity and keeps one representative per cluster, used to remove redundancy.
- **Conservation / conserved** — how unchanged a position is across related proteins; high conservation suggests functional importance.
- **Conservation grade (1–9)** — ConSurf's per-position output, from 1 (most variable) to 9 (most conserved).
- **Core peptide** — the part of the precursor that becomes the final product and carries the residues that get cross-linked.
- **DxR motif** — a conserved Asp–(any residue)–Arg pattern in graspetide cyclases; the arginine helps recognize the ring-forming acidic residue.
- **E-value** — a search statistic; the expected number of matches this good that would occur by chance. Lower means more confident.
- **Empirical Bayesian** — the statistical approach ConSurf uses to estimate how fast each position evolves; more accurate than simpler methods on small datasets.
- **Graspetide** — the modern name for the peptide family that includes microviridins, defined by ester/amide side-chain cross-links made by ATP-grasp enzymes.
- **HMM / profile HMM / HMMER / jackhmmer** — a statistical model of a protein family (HMM) and the tools that search with it; more sensitive to distant relatives than plain BLAST.
- **Homolog** — a protein related to your query by shared ancestry.
- **Leader peptide** — the part of the precursor that is recognized by the enzymes but discarded from the final product.
- **Macrocyclase / cyclase** — an enzyme that forms a ring (macrocycle) in a peptide; here, MdnB and MdnC.
- **Macrolactone (ω-ester) / macrolactam (ω-amide)** — the two cross-link types: an ester ring and an amide ring, respectively.
- **MdnA / MdnB / MdnC** — the microviridin precursor peptide (A), the lactam-forming cyclase (B), and the lactone-forming cyclase (C). Also seen as MvdE/MvdC/MvdD in other clusters.
- **MSA (multiple sequence alignment)** — related sequences arranged so equivalent positions share a column.
- **PDB / PDB ID** — the Protein Data Bank and its four-character structure codes (e.g., 5IG8).
- **Precursor peptide** — the initial gene-encoded peptide (leader + core) before modification.
- **Rate4Site** — the algorithm inside ConSurf that computes each position's evolutionary rate.
- **Residue** — a single amino acid at a specific position in a protein.
- **RiPP** — ribosomally synthesized and post-translationally modified peptide; a natural-product class made by the ribosome then chemically tailored.
- **RODEO** — a genome-mining tool that finds RiPP gene clusters using HMMs plus machine learning.
- **Serine protease / protease inhibitor** — a class of cutting enzymes, and molecules (like microviridins) that block them.
- **Sequence identity** — the percentage of positions at which two aligned sequences share the same residue.
- **SSN (sequence similarity network)** — a graph that clusters proteins by pairwise similarity to reveal family structure.
- **Twilight zone** — the ~25–35% identity range below which alignments become unreliable.
- **UniRef90** — a protein database pre-clustered so no two entries exceed 90% identity.

## References

1. Li K, Condurso HL, Li G, Ding Y, Bruner SD. Structural basis for precursor protein-directed ribosomal peptide macrocyclization. *Nature Chemical Biology*. 2016;12(11):973–979. doi:10.1038/nchembio.2200. [PDB: 5IG8, 5IG9]
2. Tietz JI, Schwalen CJ, Patel PS, Maxson T, Blair PM, Tai H-C, Zakai UI, Mitchell DA. A new genome-mining tool redefines the lasso peptide biosynthetic landscape. *Nature Chemical Biology*. 2017;13(5):470–478. doi:10.1038/nchembio.2319.
3. Lee H, Choi M, Park J-U, Roh H, Kim S. Genome mining reveals high topological diversity of ω-ester-containing peptides and divergent evolution of ATP-grasp macrocyclases. *Journal of the American Chemical Society*. 2020;142(6):3013–3023. doi:10.1021/jacs.9b12076.
4. Ramesh S, Guo X, DiCaprio AJ, De Lio AM, Harris LA, Kille BL, Pogorelov TV, Mitchell DA. Bioinformatics-guided expansion and discovery of graspetides. *ACS Chemical Biology*. 2021;16(12):2787–2797. doi:10.1021/acschembio.1c00672.
5. Walker MC, Eslami SM, Hetrick KJ, Ackenhusen SE, Mitchell DA, van der Donk WA. Precursor peptide-targeted mining of more than one hundred thousand genomes expands the lanthipeptide natural product family. *BMC Genomics*. 2020;21:387. doi:10.1186/s12864-020-06785-7.
6. Zhang et al. Large-scale biosynthetic analysis of human microbiomes reveals diverse protective ribosomal peptides. *Nature Communications*. 2025;16:3054. doi:10.1038/s41467-025-58280-w.
7. Song I, Kim S, et al. Molecular mechanism underlying substrate recognition of the peptide macrocyclase PsnB. *Nature Chemical Biology*. 2021;17:1123–1131. doi:10.1038/s41589-021-00855-x. [PDB: 7DRM, 7DRN, 7DRO, 7DRP]
8. Zhao G, Kosek D, Liu H-B, Ohlemacher SI, Blackburne B, Nikolskaya A, Makarova KS, Sun J, Barry CE III, Koonin EV, Dyda F, Bewley CA. Structural basis for a dual function ATP-grasp ligase that installs single and bicyclic ω-ester macrocycles in a new multicore RiPP natural product. *Journal of the American Chemical Society*. 2021;143(21):8056–8068. doi:10.1021/jacs.1c02316. [PDB: 7MGV]
9. Rubin GM, Patel KP, Jiang Y, Ishee AC, Seabra G, Bruner SD, Ding Y. Characterization of a dual function peptide cyclase in graspetide biosynthesis. *ACS Chemical Biology*. 2024;19(12):2525–2534. doi:10.1021/acschembio.4c00626. [PDB: 9BOU]
10. Li Y, Wang J, Zhang Z, Zhang Y, Müller R, Huo L. Deep-sea genome mining reveals cooperative ATP-grasp ligase-directed biosynthesis of pentacyclic myxomiditides with potent protease inhibition. *JACS Au*. 2026;6:607–620. doi:10.1021/jacsau.5c01626.
11. Zhang Y, et al. A distributive peptide cyclase processes multiple microviridin core peptides within a single polypeptide substrate. *Nature Communications*. 2018;9:1780.
12. Ashkenazy H, Abadi S, Martz E, Chay O, Mayrose I, Pupko T, Ben-Tal N. ConSurf 2016: an improved methodology to estimate and visualize evolutionary conservation in macromolecules. *Nucleic Acids Research*. 2016;44(W1):W344–W350. doi:10.1093/nar/gkw408.
13. Ashkenazy H, Erez E, Martz E, Pupko T, Ben-Tal N. ConSurf 2010: calculating evolutionary conservation in sequence and structure of proteins and nucleic acids. *Nucleic Acids Research*. 2010;38(Web Server):W529–W533. doi:10.1093/nar/gkq399.
14. Ben Chorin A, Masrati G, Kessel A, Narunsky A, Sprinzak J, Lahav S, Ashkenazy H, Ben-Tal N. ConSurf-DB: an accessible repository for the evolutionary conservation patterns of the majority of PDB proteins. *Protein Science*. 2020;29(1):258–267. doi:10.1002/pro.3779.
15. Using ConSurf to detect functionally important regions (ConSurf usage guidelines). *Current Protocols*. 2021. doi:10.1002/cpz1.270.
16. Mayrose I, Graur D, Ben-Tal N, Pupko T. Comparison of site-specific rate-inference methods for protein sequences: empirical Bayesian methods are superior. *Molecular Biology and Evolution*. 2004;21(9):1781–1791. doi:10.1093/molbev/msh194.
17. Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics*. 2012;28(23):3150–3152. doi:10.1093/bioinformatics/bts565.
18. Rost B. Twilight zone of protein sequence alignments. *Protein Engineering*. 1999;12(2):85–94. doi:10.1093/protein/12.2.85.
19. Hudson GA, Ramesh S. Merging RODEO with genomic enzymology workflows to guide RiPP discovery. *Methods in Enzymology*. 2025. doi:10.1016/bs.mie.2025.01.057.
20. Choi B, Link AJ. Discovery, function, and engineering of graspetides. *Trends in Chemistry*. 2023;5(8):620–633.
