# Bioinformatics Project — Master Document

*Single consolidated reference for the project: background, methodology, results so far, decisions, and next steps. This is a living document — update it as the project progresses. It is also intended as the primary source to upload to NotebookLM so the project's knowledge is stored in one place.*

*Last updated: after MdnC and MdnB ConSurf broad runs; MdnC BLAST-window curation notebook built and about to be run.*

---

## 0. Purpose & how to use this document

This is the project's memory. It captures what the project is, the biology it rests on, the computational methods and their justified parameters, the results obtained so far, the decisions made (and why), and what comes next. Detailed standalone companion files also exist (full literature review with citations; ESM methodology briefing; the analysis notebooks) — this document summarizes and links their key content so the whole project can be understood from one place.

---

## 1. Project overview

**Goal.** Build a reproducible computational pipeline that identifies functionally important residues in the microviridin/graspetide ATP-grasp macrocyclases **MdnB** and **MdnC** (and, later and separately, the precursor **MdnA**), validated against experimentally known residues, to prioritize residues worth testing by mutagenesis. The pipeline combines three complementary signals:

1. **Evolutionary conservation** (ConSurf / Rate4Site) — which positions evolution has kept unchanged.
2. **Protein language model scoring** (ESM-2, planned) — an alignment-free estimate of per-residue constraint and per-variant effect.
3. **Structure** — mapping the above onto the crystal structures and distinguishing catalytic vs recognition sites.

**Context.** This is an independent computational project within the **Bruner Lab** (University of Florida), run largely solo with mentorship from a PhD student (Amit). It wraps around the lab's experimental microviridin work, so results can feed back into the lab and serve as a strong research/portfolio artifact. The lab's own system of interest is a **multicore microviridin variant** processed by a **fused MdnBC cyclase**, with a mutant series (see §2).

**Working mode.** The researcher steers; the assistant acts as compass — giving direction, reasoning, and scaffolds, while the researcher writes the analysis and owns the work. Skills used so far: Python (Biopython, pandas basics), NCBI/BLAST, MAFFT, ChimeraX, ConSurf; GitHub for version control.

---

## 2. The biological system (background)

**RiPPs → graspetides → microviridins.** Ribosomally synthesized and post-translationally modified peptides (RiPPs) are gene-encoded peptides remodeled by tailoring enzymes. **Graspetides** (formerly "ω-ester-containing peptides") are RiPPs defined by ester and/or amide side-chain–side-chain crosslinks installed by **ATP-grasp ligases**. **Microviridins** are the founding (Group 1) graspetides — cage-like tricyclic peptides that are potent serine-protease inhibitors.

**The gene cluster and the two enzymes.** In the microviridin system the precursor **MdnA** carries an N-terminal **leader** (with the strictly conserved **PFFARFL** motif) and a C-terminal **core** (motif ~T-x-K-x-P-S-D-x-D/E-E/D) that becomes crosslinked. Two homologous ATP-grasp ligases act in order:
- **MdnC** — installs the two **ester** crosslinks (ω-esters / macrolactones; Ser/Thr donor + Asp/Glu acceptor). Also called **MvdD**.
- **MdnB** — installs the one **amide** crosslink (ω-amide / macrolactam; Lys donor + Asp/Glu acceptor). Also called **MvdC**.

**Mechanism.** Both use ATP + Mg²⁺: the enzyme phosphorylates the acidic **acceptor** (Asp/Glu) to form a reactive **acylphosphate**, which the **donor** nucleophile (Ser/Thr hydroxyl or Lys amine) attacks to form the crosslink, releasing ADP. The ATP-handling machinery is shared across the entire ATP-grasp superfamily (e.g., D-Ala–D-Ala ligase) — so those residues are *not* microviridin-specific.

**Structures (Li et al. 2016; the lab's own founding paper).** PDB **5IG9** = MdnC bound to the MdnA leader; PDB **5IG8** = MdnB (apo). Key findings:
- **Leader recognition (MdnC):** the PFFARFL leader folds into an α-helix; **MdnA Arg17** contacts **MdnC Glu191 / Asp192 / Asn195** on helix **α7** (plus Ser20–Val182). Charge reversal (E191K/D192K) abolishes binding and activity. Note: **Glu191 is often Asp in homologs** — charge conserved, identity not.
- **Conformational gate:** a **β9β10 hairpin** shifts ~25 Å between the leader-bound "open" MdnC and the "closed" MdnB, opening/closing the active site.
- **ATP/catalytic residues (MdnC), assigned by structural homology** (no nucleotide was captured in the crystals): Lys125, Lys166, Gln207, Glu215, Asp281, Glu294, Asn296.

**Comparative structures across the family** (useful reference points and validation anchors):
- **PsnB** (plesiocin; Song et al. 2021) — a conserved active-site arginine **Arg213** recognizes the ring-forming acidic acceptor.
- **PruB** (prunipeptin; Rubin et al. 2024, **from the Bruner/Ding labs at UF**) — established a conserved **DxR motif** (Asp-x-Arg; the arginine is PruB Arg233) critical for catalysis. This is the family-wide analog of PsnB Arg213. **In MdnC this is D241-W242-R243; in MdnB it is D245-W246-R247.**
- **CdnC** (chryseoviridin; Zhao et al. 2021) — dual-function ligase, leader + core + ADP complex.
- **Cooperative/multicore systems** — a distributive cyclase processes multiple cores on one precursor (Zhang et al. 2018); a deep-sea system (myxomiditides; Li et al. 2026) uses four cooperative ATP-grasp ligases, and notably some (e.g., MyxD2) **lack the E191/D192 leader-recognition residues** — direct evidence that leader recognition is lineage-specific, not universal.

**The lab's specific system (the applied target).**
- The lab studies a **multicore** microviridin precursor (many tandem K-Y-P-S-D core repeats after the PFFARFL leader) and a **fused MdnBC cyclase** (~640 aa; MdnB and MdnC domains in one chain).
- Mutant series on the fused cyclase: **H573E** (His→Glu in the MdnC-domain active site — modifies peptide first), **H193E** (His→Glu in the MdnB-domain active site), and **His→Ala** variants (not yet assayed). Hypothesis (based on side-chain pKa/charge): the Glu variant may not require ATP; the His (wild type) requires some ATP; the Ala variant is inactive regardless of ATP. ITC assays with wild-type MdnBC + precursor + ATP have been run; the mutants are the experimental frontier and the eventual ground truth for calibrating predictions.
- **Lab-provided sequences (stored here so they never need re-pasting):**

  *Fused MdnBC cyclase (~640 aa; MdnB domain N-terminal, MdnC domain C-terminal):*
  ```
  MSVSLVLLLTHRRDSFTIDRVHAAVERLGATAVRVDTDRFPAELSLGVRITGGLLTGRLHGDGRTIDLSEVGAVWMRRLWPPGGLEALAPRWRGTSHHHSHLALTQLLPLLSHARWLDRIDRQLAAESKPRQLVEAAQVGLTVPDTLITNASGEVQRFDREQGPLVTKLLEPIAYGMEGGRGDFMYTSRVTTHDLAAMDGLRWVPQIFQPEVPKAKELRVVVVGEQVFAGAIDTGASTRGTVDWRRLRANEGPPWEPAQLPPPVVTATRALLARFDLTFGVLDFIVTPAGEYVFLDLNPAGEWGWLERDLDLPISEAIARWLVDAASAATNDEDPERSASFAEASTDLHAVPRPPAVLPSTSTPAAGPASADHHGPPAPAVLRDPTTRPAVLIITHSGDNECIDTVSTAIRNRGGRPIRLDTDRYPTSVGLSTLQGQSSTGLVIAGDEPVPLESLQAVWYRRFAAGGRLPDSMGDTREAAVGETRRTLQGTIANLGCFQLDPLAAVQRTHHKELQLRLAAEEGLPVPRTLVTNRAQDAREFWRRLDGRVVTKMQHSFAIYREGRETVVFTSRVREHDLDALDGLRLCPMTFQEEIPKALELRVTIVGTRVMTAAIDSARRAKTEVDWRRDGRGLLHDWDPYELPPPVARRMLALQRRLGLNYGAADMILTPDGRHVLLEVNPVGEFFWLDRDPGLPISDSIAAVLLGQAERNVDRLGG
  ```

  *Multicore MdnA precursor (PFFARFL leader + tandem K-Y-P-S-D core repeats):*
  ```
  MAAQLMNEVAPALLKEQPFFARFLEDQEDQSDDEDPADAGGTPTTLKYPSDEEDGLAQTKKYPSDHEDGLVQTRKYPSDEEDSPVQTKKYPSDHEDGPVHTMKYPSDEEDGPVQTKKYPSDHEDGIVQTRKYPSDHEDGIAQTKKYPSDEEDSAAVTEKFPSDHEDGIAYTMKFPSDWEDGTPGR
  ```
- A parallel lab project engineers the pathway onto **yeast surface display** (Aga1/Aga2 + SpyTag/SpyCatcher, co-displaying MdnA with MdnBC) toward drug-targeting applications — context, not the focus of this computational project.

---

## 3. Validation residue set (the "answer key")

These are the experimentally supported or family-conserved residues a correct conservation/ESM map should recover. They are the benchmark: if the pipeline flags these, its predictions on *untested* residues can be trusted.

| Residue(s) | Role | Source |
|---|---|---|
| MdnC K125, K166, Q207, E215, D281, E294, N296 | ATP/metal binding (catalytic core) | Li et al. 2016 (structural assignment) |
| MdnC E191, D192, N195 (α7) | Leader-peptide recognition | Li et al. 2016 (E191K/D192K abolishes binding/activity) |
| DxR arginine — MdnC R243 (D241xR243); MdnB R247 (D245xR247) | Acceptor recognition / catalysis | Rubin et al. 2024 (PruB); Song et al. 2021 (PsnB Arg213) |
| Fused MdnBC H573E (MdnC domain), H193E (MdnB domain), His→Ala | Active-site His; ATP-dependence hypothesis | Lab's own experiments |

---

## 4. Methodology

### 4a. Conservation pipeline (ConSurf) — parameters and rationale

Run **MdnB and MdnC separately** — they are divergent clades (a lactam-former vs a lactone-former); pooling them averages away the very residues that distinguish them.

| Step | Choice | Rationale (source) |
|---|---|---|
| Homolog collection | HMMER/UniRef90 (ConSurf default) **or** curated BLAST set | Diversity matters more than raw count |
| Identity window | ~**35–90%** to query | <35% = "twilight zone," unreliable alignment (Rost 1999); >90–95% = redundant |
| Redundancy removal | **CD-HIT 0.9** (or equivalent) | Removes over-sequenced-organism / near-duplicate bias (Fu et al. 2012) |
| Homolog count | **50–300**, diverse | ConSurf needs ≥50 for stable rates; >300 diminishing returns |
| Alignment | **MAFFT** | ConSurf default; high accuracy |
| Scoring | ConSurf / **Rate4Site**, empirical Bayesian → grades **1–9** (9 = most conserved) | Mayrose et al. 2004; Ashkenazy et al. 2016 |
| Interpretation | Map onto 5IG8/5IG9; prioritize conserved residues **outside** the universal ATP pocket (leader/substrate clefts) | ATP-pocket residues are superfamily-wide and uninformative for family specificity |

### 4b. ESM-2 protein language model (planned layer)

**Model/scheme:** ESM-2 **650M**, **masked-marginal** scoring as primary; cross-check specific mutants with the ESM-1v five-model ensemble. (8M/35M run on CPU/M1 for quick tests; 650M is the accuracy/cost sweet spot, best on a free Colab GPU.)

**Two products:** (1) a **per-residue constraint track** (aggregate the log-likelihood ratios over all 19 substitutions, or use predictive entropy) — directly comparable to ConSurf grades; (2) **per-variant effect scores** for specific mutations (validate against the residue set in §3).

**Validation logic:** functional residues should score as strongly constrained. Where ESM and ConSurf **agree** = high confidence. Where they **disagree** is informative: conserved-but-ESM-low → likely a family-specific catalytic residue conservation catches; ESM-constrained-but-ConSurf-low → possible shallow-alignment artifact, ESM adds info.

**Later:** add a structure-aware track (SaProt or ESM-IF1 on 5IG9/5IG8; AlphaFold-Multimer/AF3 for enzyme–precursor complexes of homologs), and once ~20–50 mutant activities exist, fuse tracks with a simple Hsu-style augmented ridge regression rather than a bespoke model.

**Known pitfalls (state these in any writeup):**
- ESM scores **constraint/plausibility, not mechanism** — it will likely *not* predict that H573E removes the ATP requirement; that is a mechanistic effect neither conservation nor a language model is built to capture.
- The **fused MdnBC** construct is partly out-of-distribution — score full-length *and* per-domain, and compare to native standalone MdnB/MdnC.
- The **multicore MdnA** precursor is short and repetitive — score the whole precursor and individual leader+core units separately.
- RiPP-specific PLM precedent is thin (LassoESM, LazBF/LazDEF on other RiPP classes); there is essentially no prior PLM variant-effect study on graspetide macrocyclases — so this application is novel, and transfer rests on the general variant-effect literature.

---

## 5. Results so far

### 5a. MdnC ConSurf — broad run

**Flagged-residue grades (1–9):**
- ATP/catalytic core — **all grade 9**: K125, K166, Q207, E215, D281, E294, N296. (E215 is grade 9 but the column is ~60% D / 38% E — conserved *as an acidic residue*, not strictly E.)
- **DxR arginine R243 — grade 9** (D241 grade 9, W242 grade 8).
- α7 leader cluster — **low**: **E191 = grade 3, D192 = grade 4, N195 = grade 5** (N195 had ~half-coverage in the MSA).

**Homolog set diagnosis (the key finding):** the set was **too broad** — 300 sequences, **median only ~28% identity** to MdnC, with **280/300 below 35%**. Headers were dominated by distant Actinobacteria ATP-grasp ligases (*Streptomyces, Frankia, Nonomuraea, Sphaerisporangium, Tenggerimyces, Anaeromyxobacter*); only the top hit was *Microcystis* (a genuine microviridin producer). ConSurf's own 35%-identity setting did not exclude these because it appears to measure identity over the HMMER-aligned core region, where distant enzymes still clear 35%.

**Interpretation:** the result is coherent and validates the pipeline — the universal catalytic core and the DxR arginine conserve strongly because they are shared across the whole ATP-grasp/graspetide superfamily. The α7 leader cluster washes out because leader recognition is **microviridin-specific**, and ~93% of these homologs are not microviridins (E191's column even shows 17% Ala / 9% Gly — non-charged residues no leader-recognizer would tolerate). So the broad set can only ever reveal superfamily-wide residues; it cannot surface the microviridin-specific residues that are the project's actual target.

### 5b. MdnB ConSurf — broad run

- **DxR arginine R247 — grade 9** (D245 grade 9, W246 grade 8); catalytic core conserved. Same clean validation as MdnC.
- Same **too-broad** homolog set (median ~29% identity; ~275/299 below 35%).
- **Numbering caution:** MdnB positions 191–195 are R/V/K/A/E — **not** the leader-recognition cluster (that is MdnC's numbering). MdnB's leader-recognition equivalent must be found by **structural superposition** of 5IG8 onto 5IG9, not by matching residue numbers. MdnB also binds the leader ~10× more weakly than MdnC, so its recognition residues are expected to be less conserved regardless.

### 5c. The two-tier framing (turns the problem into a design)

- The **broad run** is a valid **superfamily-wide map**: it authoritatively confirms the universal catalytic core + DxR. Keep it.
- A **focused microviridin run** is needed to surface the **family-specific** layer (leader/substrate recognition).
- Comparing the two directly shows which residues are universal vs microviridin-specific — itself a useful analytical result, and exactly what is needed to prioritize mutagenesis.

---

## 6. Decisions log & current state

**Decisions made (and why):**
- **Analyze MdnB and MdnC separately** — divergent clades; pooling averages away distinguishing residues.
- **Do not use the lab's private anti-SMASH deep-sea database** — unsorted, function-undetermined, biased, and very large; it would inject unrelated sequences with no way to filter by function. Public databases are cleaner and defensible. (Keep it in mind as a *future discovery* project, not an input here.)
- **Precursor (MdnA) is not analyzed with ConSurf** — it's a short, hypervariable, multicore peptide; ConSurf is ill-suited. Handle later via targeted alignment of known graspetide precursors and/or the ESM track.
- **Fix the too-broad set by curating a microviridin-focused homolog set**, not by nudging ConSurf's identity slider (which measures identity over the aligned core and can't cleanly separate microviridin macrocyclases from distant ATP-grasp cousins).
- **Curation route:** start with a **BLAST-window** approach (fast, query-centered), then optionally **enrich with the RODEO graspetide dataset** (function-verified). Combining both is the strongest option; BLAST-window alone is a legitimate first pass.
- **MdnC first**, because the whole validation set lives there.

**Current task:** running the **MdnC BLAST-window curation notebook** (`MdnC_BLASTwindow_curation.ipynb`) to build a curated 35–90% microviridin-macrocyclase MSA, then re-running ConSurf on that custom MSA to test whether the α7 cluster rises off the floor.

**Immediate success test:** in the focused run, does **E191/D192/N195** climb toward conserved (grade ≥6, expected as acidic E/D — not necessarily 9), while the catalytic core and DxR R243 stay high? If yes, the focused set fixed the signal and new microviridin-specific positions can be read off. Watch the final sequence count — if <~50, real microviridin macrocyclases in-band are scarce, and the next move is RODEO enrichment.

---

## 7. Homolog curation approach (detail)

**BLAST-window (first pass):**
1. BLAST the query (MdnC, then MdnB) against NCBI nr.
2. Keep hits **35–90% identical to the query**, measured as *fraction of the full query that matches* (stricter than BLAST local identity — this is what prevents distant enzymes from leaking in), with ≥70% coverage.
3. **Macrocyclase screen:** keep ~280–360 aa sequences (full-domain ATP-grasp ligases, not fragments/fusions).
4. **Redundancy removal:** drop sequences >90% identical to each other (CD-HIT or equivalent).
5. Target 50–150 diverse sequences; MAFFT-align; upload to ConSurf as a custom MSA (with 5IG9 chain A / 5IG8 chain A as the structure).

**RODEO enrichment (strengthen / fall back):** download the published RODEO graspetide dataset (Ramesh et al. 2021 — ~3,900 clusters, pre-sorted into 24 groups), pull the **Group 1 (microviridin) macrocyclases**, split into lactone-formers (MdnC/MvdD-type) and lactam-formers (MdnB/MvdC-type). This gives functional certainty and clean group labels; it is a fixed 2021 snapshot and must be mapped to the identity window relative to the query.

**Best result = combine:** RODEO backbone (function-verified, correctly typed) + BLAST-window hits (query-centered, includes post-2021 sequences), merged, de-duplicated, CD-HIT 0.9, in the 35–90% band.

---

## 8. Files & artifacts inventory

- **Literature review** — full background with verified citations (graspetide biology, conservation methodology). *(standalone file)*
- **ESM-2 methodology briefing** — protein language model methods, RiPP-specific PLM landscape, concrete plan for MdnB/MdnC. *(standalone file)*
- **This master document** — consolidated project reference.
- **Notebooks:** Week-2 homolog + alignment notebooks (superseded); **`MdnC_BLASTwindow_curation.ipynb`** (current — curates the focused MdnC set).
- **ConSurf outputs:** MdnC broad-run grades + MSA; MdnB broad-run grades + MSA.
- **Structures:** PDB 5IG9 (MdnC + MdnA leader), 5IG8 (MdnB), local copies in `data/`.
- **Lab-provided sequences:** the fused MdnBC cyclase and the multicore MdnA precursor.
- **GitHub repo** — the working project folder (notebooks/, data/, figures/, results/), synced via GitHub Desktop.

---

## 9. Open questions & next steps

1. **Does the focused MdnC set fix α7?** (immediate test — running now).
2. **Are there ≥50 microviridin macrocyclases in the 35–90% band?** If not, enrich with RODEO or widen slightly.
3. **MdnB leader-recognition residues** — identify by structural superposition of 5IG8 onto 5IG9 in ChimeraX (not by number).
4. **ESM-2 implementation** — per-residue and per-variant tracks; validate against §3; compare with ConSurf.
5. **Precursor (MdnA)** — separate sub-analysis via targeted graspetide-precursor alignment and/or ESM (handle the multicore repeats carefully).
6. **Fused MdnBC handling** — score full-length vs per-domain; the H573E/H193E mutants are the eventual mechanistic ground truth.
7. **Two-tier comparison** — formally contrast the broad (superfamily) and focused (microviridin) maps to classify residues as universal vs family-specific.

---

## 10. References (verified)

1. Li K, Condurso HL, Li G, Ding Y, Bruner SD. Structural basis for precursor protein-directed ribosomal peptide macrocyclization. *Nat Chem Biol.* 2016;12(11):973–979. doi:10.1038/nchembio.2200. [PDB 5IG8, 5IG9]
2. Tietz JI, et al. A new genome-mining tool redefines the lasso peptide biosynthetic landscape. *Nat Chem Biol.* 2017;13(5):470–478. doi:10.1038/nchembio.2319.
3. Lee H, Choi M, Park J-U, Roh H, Kim S. Genome mining reveals high topological diversity of ω-ester-containing peptides and divergent evolution of ATP-grasp macrocyclases. *J Am Chem Soc.* 2020;142(6):3013–3023. doi:10.1021/jacs.9b12076.
4. Ramesh S, et al. Bioinformatics-guided expansion and discovery of graspetides. *ACS Chem Biol.* 2021;16(12):2787–2797. doi:10.1021/acschembio.1c00672.
5. Walker MC, et al. Precursor peptide-targeted mining of more than one hundred thousand genomes expands the lanthipeptide natural product family. *BMC Genomics.* 2020;21:387. doi:10.1186/s12864-020-06785-7.
6. Song I, et al. Molecular mechanism underlying substrate recognition of the peptide macrocyclase PsnB. *Nat Chem Biol.* 2021;17:1123–1131. doi:10.1038/s41589-021-00855-x. [PDB 7DRM/7DRN/7DRO/7DRP]
7. Zhao G, et al. Structural basis for a dual function ATP-grasp ligase that installs single and bicyclic ω-ester macrocycles in a new multicore RiPP natural product. *J Am Chem Soc.* 2021;143(21):8056–8068. doi:10.1021/jacs.1c02316. [PDB 7MGV]
8. Rubin GM, Patel KP, Jiang Y, Ishee AC, Seabra G, Bruner SD, Ding Y. Characterization of a dual function peptide cyclase in graspetide biosynthesis. *ACS Chem Biol.* 2024;19(12):2525–2534. doi:10.1021/acschembio.4c00626. [PDB 9BOU]
9. Li Y, et al. Deep-sea genome mining reveals cooperative ATP-grasp ligase-directed biosynthesis of pentacyclic myxomiditides. *JACS Au.* 2026;6:607–620. doi:10.1021/jacsau.5c01626.
10. Zhang Y, et al. A distributive peptide cyclase processes multiple microviridin core peptides within a single polypeptide substrate. *Nat Commun.* 2018;9:1780.
11. Choi B, Link AJ. Discovery, function, and engineering of graspetides. *Trends Chem.* 2023;5(8):620–633.
12. Ashkenazy H, et al. ConSurf 2016. *Nucleic Acids Res.* 2016;44(W1):W344–W350. doi:10.1093/nar/gkw408.
13. Ben Chorin A, et al. ConSurf-DB. *Protein Sci.* 2020;29(1):258–267. doi:10.1002/pro.3779.
14. Mayrose I, Graur D, Ben-Tal N, Pupko T. Comparison of site-specific rate-inference methods; empirical Bayesian methods are superior. *Mol Biol Evol.* 2004;21(9):1781–1791. doi:10.1093/molbev/msh194.
15. Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics.* 2012;28(23):3150–3152. doi:10.1093/bioinformatics/bts565.
16. Rost B. Twilight zone of protein sequence alignments. *Protein Eng.* 1999;12(2):85–94. doi:10.1093/protein/12.2.85.
17. Meier J, et al. Language models enable zero-shot prediction of the effects of mutations on protein function. *NeurIPS.* 2021. doi:10.1101/2021.07.09.450648. [ESM-1v]
18. Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science.* 2023;379:1123–1130. doi:10.1126/science.ade2574. [ESM-2]
19. Notin P, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *NeurIPS Datasets & Benchmarks.* 2023.
20. Hsu C, et al. Learning protein fitness models from evolutionary and assay-labeled data. *Nat Biotechnol.* 2022;40:1114–1122. doi:10.1038/s41587-021-01146-5.
21. Mi X, et al. LassoESM: a domain-adapted protein language model for lasso peptides. *Nat Commun.* 2025. doi:10.1038/s41467-025-63412-3.

*(Full annotated citations and a glossary are in the standalone literature-review and ESM-methodology files.)*
