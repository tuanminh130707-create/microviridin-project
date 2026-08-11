# Bioinformatics Project — Master Document

*Single consolidated reference for the project: background, methodology, results so far, decisions, and next steps. This is a living document — update it as the project progresses. It is also intended as the primary source to upload to NotebookLM so the project's knowledge is stored in one place.*

*Last updated: **MdnC focused pipeline complete through ConSurf.** A curated, decontaminated, microviridin-specific MdnC homolog set was built (HMM classification of RODEO → taxonomic filter → BLAST-window merge → CD-HIT → IQ-TREE clade), and the focused ConSurf run is done. Headline result: **Asp192 is an invariant, ConSurf-"functional" grade-9 residue** (it rose from grade 4 in the broad run), while Glu191 and Asn195 are graded lower for reasons that are biologically informative. See §5d–§5f.*

---

## 0. Purpose & how to use this document

This is the project's memory. It captures what the project is, the biology it rests on, the computational methods and their justified parameters, the results obtained so far, the decisions made (and why), and what comes next. Detailed standalone companion files also exist (full literature review with citations; ESM methodology briefing; the **ConSurf homolog-curation methodology review**; the analysis notebooks) — this document summarizes and links their key content so the whole project can be understood from one place.

**Orientation for a returning reader:** the project passed a major turning point. The *broad* ConSurf runs (§5a–5c) correctly recovered the universal catalytic core but washed out the microviridin-specific leader-recognition residues — the actual target. The fix was a full homolog-curation pipeline (§4c, §5d) that turned out to require more than the originally planned identity window: it also required **HMM family classification and a taxonomic filter**, because identity alone cannot separate microviridin macrocyclases from the many non-cyanobacterial ester-forming graspetide ligases that share the fold. That pipeline is now complete for MdnC and produced a clean, decontaminated 71-sequence set, and the focused ConSurf on it (§5e) delivered the primary result.

---

## 1. Project overview

**Goal.** Build a reproducible computational pipeline that identifies functionally important residues in the microviridin/graspetide ATP-grasp macrocyclases **MdnB** and **MdnC** (and, later and separately, the precursor **MdnA**), validated against experimentally known residues, to prioritize residues worth testing by mutagenesis. The pipeline combines three complementary signals:

1. **Evolutionary conservation** (ConSurf / Rate4Site) — which positions evolution has kept unchanged.
2. **Protein language model scoring** (ESM-2, planned) — an alignment-free estimate of per-residue constraint and per-variant effect.
3. **Structure** — mapping the above onto the crystal structures and distinguishing catalytic vs recognition sites.

**Context.** This is an independent computational project within the **Bruner Lab** (University of Florida), run largely solo with mentorship from a PhD student (Amit). It wraps around the lab's experimental microviridin work, so results can feed back into the lab and serve as a strong research/portfolio artifact. The lab's own system of interest is a **multicore microviridin variant** processed by a **fused MdnBC cyclase**, with a mutant series (see §2).

**Working mode.** The researcher steers; the assistant acts as compass — giving direction, reasoning, and scaffolds, while the researcher writes the analysis and owns the work. Skills used so far: Python (Biopython, pandas), NCBI/BLAST, **HMMER (hmmscan/hmmpress)**, **CD-HIT**, MAFFT, **trimAl**, **IQ-TREE 2 (ModelFinder, UFBoot, SH-aLRT)**, **FastTree**, ChimeraX, ConSurf; GitHub for version control.

---

## 2. The biological system (background)

**RiPPs → graspetides → microviridins.** Ribosomally synthesized and post-translationally modified peptides (RiPPs) are gene-encoded peptides remodeled by tailoring enzymes. **Graspetides** (formerly "ω-ester-containing peptides") are RiPPs defined by ester and/or amide side-chain–side-chain crosslinks installed by **ATP-grasp ligases**. **Microviridins** are the founding (Group 1) graspetides — cage-like tricyclic peptides that are potent serine-protease inhibitors.

**The gene cluster and the two enzymes.** In the microviridin system the precursor **MdnA** carries an N-terminal **leader** (with the strictly conserved **PFFARFL** motif) and a C-terminal **core** (motif ~T-x-K-x-P-S-D-x-D/E-E/D) that becomes crosslinked. Two homologous ATP-grasp ligases act in order:
- **MdnC** — installs the two **ester** crosslinks (ω-esters / macrolactones; Ser/Thr donor + Asp/Glu acceptor). Also called **MvdD**.
- **MdnB** — installs the one **amide** crosslink (ω-amide / macrolactam; Lys donor + Asp/Glu acceptor). Also called **MvdC**.

> **Naming map (used constantly downstream, and a genuine trap):** the two are cross-named. Mvd**D** = Mdn**C** = the **ester**-former; Mvd**C** = Mdn**B** = the **amide**-former. The TIGRFAM/NCBIfam HMM profiles carry these as their *accessions* — **TIGR04184 = ATPgraspMvdD (ester/MdnC-type)** and **TIGR04185 = ATPgraspMvdC (amide/MdnB-type)** — but each model's internal **NAME** is `ATPgraspMvdD` / `ATPgraspMvdC`, *not* the TIGR accession. All HMM parsing keys on the NAME.

**Mechanism.** Both use ATP + Mg²⁺: the enzyme phosphorylates the acidic **acceptor** (Asp/Glu) to form a reactive **acylphosphate**, which the **donor** nucleophile (Ser/Thr hydroxyl or Lys amine) attacks to form the crosslink, releasing ADP. The ATP-handling machinery is shared across the entire ATP-grasp superfamily (e.g., D-Ala–D-Ala ligase) — so those residues are *not* microviridin-specific.

**Structures (Li et al. 2016; the lab's own founding paper).** PDB **5IG9** = MdnC bound to the MdnA leader; PDB **5IG8** = MdnB (apo). Key findings:
- **Leader recognition (MdnC):** the PFFARFL leader folds into an α-helix; **MdnA Arg17** contacts **MdnC Glu191 / Asp192 / Asn195** on helix **α7** (plus Ser20–Val182). Charge reversal (E191K/D192K) **abolishes** binding and activity; the milder E191A/D192A **retains partial cyclization** but loses tight binding. Note: **Glu191 is often Asp in homologs** — charge conserved, identity not.
- **Conformational gate:** a **β9β10 hairpin** shifts ~25 Å between the leader-bound "open" MdnC and the "closed" MdnB, opening/closing the active site.
- **ATP/catalytic residues (MdnC), assigned by structural homology** (no nucleotide was captured in the crystals): Lys125, Lys166, Gln207, Glu215, Asp281, Glu294, Asn296.

**Comparative structures across the family** (useful reference points and validation anchors):
- **PsnB** (plesiocin; Song et al. 2021) — a conserved active-site arginine **Arg213** recognizes the ring-forming acidic acceptor.
- **PruB** (prunipeptin; Rubin et al. 2024, **from the Bruner/Ding labs at UF**) — established a conserved **DxR motif** (Asp-x-Arg; the arginine is PruB Arg233) critical for catalysis. This is the family-wide analog of PsnB Arg213. **In MdnC this is D241-W242-R243; in MdnB it is D245-W246-R247.**
- **CdnC** (chryseoviridin; Zhao et al. 2021) — dual-function ligase, leader + core + ADP complex; a *Chryseobacterium* (Bacteroidota) enzyme. **Relevant to §5d:** this and related non-cyanobacterial ester-formers are precisely the contaminants the curation pipeline had to exclude.
- **Cooperative/multicore systems** — a distributive cyclase processes multiple cores on one precursor (Zhang et al. 2018); a deep-sea system (myxomiditides; Li et al. 2026) uses four cooperative ATP-grasp ligases, and notably some (e.g., MyxD2) **lack the E191/D192 leader-recognition residues** — direct evidence that leader recognition is lineage-specific, not universal. This prediction was independently borne out by the conservation result (§5e): the leader-recognition residues are microviridin-specific.

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

| Residue(s) | Role | Source | Recovered? |
|---|---|---|---|
| MdnC K125, K166, Q207, E215, D281, E294, N296 | ATP/metal binding (catalytic core) | Li et al. 2016 (structural assignment) | ✅ broad run, all grade 9 |
| MdnC E191, D192, N195 (α7) | Leader-peptide recognition | Li et al. 2016 (E191K/D192K abolishes binding/activity) | ✅ **focused run** — see §5e (D192 grade 9; E191/N195 graded lower, and the reason is informative) |
| DxR arginine — MdnC R243 (D241xR243); MdnB R247 (D245xR247) | Acceptor recognition / catalysis | Rubin et al. 2024 (PruB); Song et al. 2021 (PsnB Arg213) | ✅ broad run, grade 9 |
| Fused MdnBC H573E (MdnC domain), H193E (MdnB domain), His→Ala | Active-site His; ATP-dependence hypothesis | Lab's own experiments | pending ESM / experimental |

**Reading the α7 recovery honestly:** the original expectation was that all three α7 residues would score conserved. The focused run supports this **strongly for Asp192** (grade 9, invariant) but shows **Glu191 and Asn195 as more tolerant** — which, rather than a failure, is a refinement that matches the mutagenesis (see §5e). The success criterion "did the focused set lift the α7 signal off the floor?" is met: **Asp192 rose from grade 4 (broad) to grade 9 (focused).**

---

## 4. Methodology

### 4a. Conservation pipeline (ConSurf) — parameters and rationale

Run **MdnB and MdnC separately** — they are divergent clades (a lactam-former vs a lactone-former); pooling them averages away the very residues that distinguish them.

| Step | Choice | Rationale (source) |
|---|---|---|
| Homolog set | **Curated, provided as an external MSA** (see §4c) — *not* ConSurf's built-in HMMER/UniRef90 search | The built-in search re-introduces exactly the distant-homolog contamination the curation removes (§5d). Uploading a custom MSA bypasses it. |
| Identity window | ~**35–90%** to query for BLAST-window collection | <35% = "twilight zone," unreliable alignment (Rost 1999); >90–95% = redundant |
| Redundancy removal | **CD-HIT 0.90** | 0.80–0.95 is the ConSurf-recommended band; below 0.80 strips diversity (Rubin 2021; Fu et al. 2012) |
| Homolog count | **50–150** (ConSurf caps at 300), diverse | ConSurf needs enough for stable rates; the final MdnC set is **71** |
| Alignment | **MAFFT L-INS-i** | Highest-accuracy MAFFT mode; feeds both tree and ConSurf |
| Scoring | ConSurf / **Rate4Site**, **empirical Bayesian** → grades **1–9** (9 = most conserved), with the "insufficient data" flag | Mayrose et al. 2004 (EB superior at small N); Ashkenazy et al. 2016 |
| Substitution model | **LG** (chosen by ModelFinder on the tree; ConSurf applies gamma rates on top) | Keeps tree and conservation on consistent evolutionary assumptions |
| Query / structure | Query = **WP_012265553** (the *Microcystis* MdnC), mapped onto **5IG9 chain A** | The query must be the *Microcystis* sequence so grades map onto the 5IG9 structure at PDB numbering |
| Interpretation | Map onto 5IG8/5IG9; prioritize conserved residues **outside** the universal ATP pocket (leader/substrate clefts) | ATP-pocket residues are superfamily-wide and uninformative for family specificity |

### 4b. ESM-2 protein language model (planned layer)

**Model/scheme:** ESM-2 **650M**, **masked-marginal** scoring as primary; cross-check specific mutants with the ESM-1v five-model ensemble. (8M/35M run on CPU/M1 for quick tests; 650M is the accuracy/cost sweet spot, best on a free Colab GPU.)

**Two products:** (1) a **per-residue constraint track** (aggregate the log-likelihood ratios over all 19 substitutions, or use predictive entropy) — directly comparable to ConSurf grades; (2) **per-variant effect scores** for specific mutations (validate against the residue set in §3).

**Validation logic:** functional residues should score as strongly constrained. Where ESM and ConSurf **agree** = high confidence. Where they **disagree** is informative: conserved-but-ESM-low → likely a family-specific catalytic residue conservation catches; ESM-constrained-but-ConSurf-low → possible shallow-alignment artifact, ESM adds info. **A concrete test now available:** does ESM-2 rank **D192 ≫ E191 ≈ N195** in mutational sensitivity, as both ConSurf (§5e) and the Li 2016 mutagenesis do? Agreement there would be a strong third independent line.

**Later:** add a structure-aware track (SaProt or ESM-IF1 on 5IG9/5IG8; AlphaFold-Multimer/AF3 for enzyme–precursor complexes of homologs), and once ~20–50 mutant activities exist, fuse tracks with a simple Hsu-style augmented ridge regression rather than a bespoke model.

**Known pitfalls (state these in any writeup):**
- ESM scores **constraint/plausibility, not mechanism** — it will likely *not* predict that H573E removes the ATP requirement; that is a mechanistic effect neither conservation nor a language model is built to capture.
- The **fused MdnBC** construct is partly out-of-distribution — score full-length *and* per-domain, and compare to native standalone MdnB/MdnC.
- The **multicore MdnA** precursor is short and repetitive — score the whole precursor and individual leader+core units separately.
- RiPP-specific PLM precedent is thin; there is essentially no prior PLM variant-effect study on graspetide macrocyclases — so this application is novel.

### 4c. Focused microviridin homolog curation — HMM classification → taxonomic filter → phylogenetic clade *(NEW — the methodology that actually surfaces the family-specific signal)*

The broad ConSurf run failed for a specific, diagnosable reason (§5a): a homolog set defined by identity alone cannot separate microviridin macrocyclases from the far larger pool of ester-forming graspetide ligases that share the ATP-grasp fold across many bacterial phyla. Fixing this took a four-stage curation pipeline. **Each stage answers a different question, and all four are necessary** — the central lesson is that no single filter suffices.

**Stage 1 — HMM family classification (separates ester- from amide-formers).**
- Combine two curated profile HMMs into one searchable database: **TIGR04184** (`ATPgraspMvdD`, ester/MdnC-type, gathering threshold **475 bits**) and **TIGR04185** (`ATPgraspMvdC`, amide/MdnB-type, GA **425 bits**), then `hmmpress`.
- Run `hmmscan --cut_ga` (use the curators' gathering thresholds, not a generic E-value) and classify each sequence by which model it clears: MvdD-only → **MdnC-type**; MvdC-only → **MdnB-type**; both → ambiguous (excluded); neither → not a member.
- *Why HMM and not pairwise identity:* pairwise BLOSUM62 identity gives a **non-bimodal** distribution for this family (the discriminating features are locally concentrated), so identity can't cleanly split the two chemistries; the split HMMs can. Validated on positive controls (WP_012265553 ester → MvdD by ~272 bits; WP_012265554 amide → MvdC by ~226 bits).

**Stage 2 — taxonomic filter to Cyanobacteria (separates microviridin from other-phylum ester-formers).**
- **This is the stage the original plan lacked, and the key methodological finding of the project.** Clearing the ester-former HMM threshold is *necessary but not sufficient*: it certifies "this is an ester-forming graspetide ligase," **not** "this is a microviridin." Microviridin is a cyanobacterial natural product, but ester-forming graspetide ligases occur across Bacteroidota, Pseudomonadota (Proteobacteria), Acidobacteriota, and Myxococcota — and these pass the MvdD HMM. They must be removed by **taxonomy**, using the phylum annotation in the RODEO spreadsheet (keep `Phylum == Cyanobacteria`).
- *Why taxonomy and not a tighter HMM or identity cutoff:* HMM gathering thresholds classify by fold/family, not by leader-recognition specificity; and the non-cyanobacterial ester-formers sit on **both sides** of the identity/bit-score line, so no threshold cleanly separates them (§5d).

**Stage 3 — supplement with a query-centered BLAST-window set (adds post-2021 diversity), then merge + deduplicate.**
- RODEO is a fixed 2021 snapshot; a focused BLAST of MdnC vs current `nr`, filtered to the 35–90% window and passed through the *same* HMM + taxonomy gates, adds sequences RODEO cannot contain. Merge with the cyanobacterial RODEO set, dedup on accession, then **CD-HIT 0.90** to collapse near-duplicates. **Re-insert the query anchor after CD-HIT** (clustering can discard it).

**Stage 4 — phylogenetic clade definition + readiness check (the final set-definition step).**
- Align the set + a proximal non-cyanobacterial outgroup (MAFFT L-INS-i), **trim a copy for the tree only** (trimAl `-gt 0.2`; never trim the ConSurf MSA — it would delete the target columns), build the tree (FastTree preview → **IQ-TREE 2** with ModelFinder `-m MFP`, `-B 1000` UFBoot, `-alrt 1000` SH-aLRT), and confirm the cyanobacterial in-group holds together with the anchor central.
- The ConSurf MSA is the **untrimmed alignment restricted to the clade, outgroup removed, anchor included**. A **readiness diagnostic** confirms before submission: set size in band, anchor present and full-length, **≥90% occupancy at the target columns 191/192/195** (a gappy target column causes a false "insufficient data" flag), and no non-cyanobacterial leaks.

*(The parameter choices above — 50–150 sequences, CD-HIT 0.80–0.95, empirical-Bayesian rates, external-MSA upload, monophyletic-clade set definition, vetting the outgroup for compositional bias — are grounded in a standalone **ConSurf homolog-curation methodology review** compiled for this project, which also confirms the graded α7 pattern against Li 2016.)*

---

## 5. Results so far

### 5a. MdnC ConSurf — broad run *(historical; the run that motivated the whole curation effort)*

**Flagged-residue grades (1–9):**
- ATP/catalytic core — **all grade 9**: K125, K166, Q207, E215, D281, E294, N296. (E215 is grade 9 but ~60% D / 38% E — conserved *as an acidic residue*.)
- **DxR arginine R243 — grade 9** (D241 grade 9, W242 grade 8).
- α7 leader cluster — **low**: **E191 = grade 3, D192 = grade 4, N195 = grade 5** (N195 had ~half-coverage in the MSA).

**Homolog set diagnosis:** too broad — 300 sequences, **median only ~28% identity** to MdnC, **280/300 below 35%**. The taxonomy tells the whole story in one ratio: **145/300 sequences were *Streptomyces* (~48%)** and the pool was overwhelmingly Actinobacterial (Kitasatospora 16, Actinomadura 9, Nonomuraea 6, plus Nocardia, Micromonospora, Kribbella, Streptosporangium…), with only **~5 cyanobacterial sequences (1 *Microcystis* + 4 *Nostoc*)** — i.e. essentially **one genuine microviridin producer in the entire set** (the top hit, *Microcystis* B0JH80). ConSurf's 35%-identity setting didn't exclude the distant enzymes (it measures identity over the HMMER-aligned core, where they still clear 35%). *(Verified directly from the recovered broad-run MSA, `msa_fasta.aln`.)*

**Interpretation:** coherent and pipeline-validating — the universal catalytic core and DxR conserve because they are shared across the whole ATP-grasp superfamily; the α7 cluster washes out because leader recognition is microviridin-specific and ~98% of the set are non-cyanobacterial (only one is a known microviridin producer). The broad set can only reveal superfamily-wide residues.

### 5b. MdnB ConSurf — broad run *(historical)*

- **DxR arginine R247 — grade 9** (D245 grade 9, W246 grade 8); catalytic core conserved. Same clean validation as MdnC.
- Same **too-broad** homolog set (median ~29% identity; ~275/299 below 35%).
- **Numbering caution:** MdnB positions 191–195 are R/V/K/A/E — **not** the leader-recognition cluster (that is MdnC's numbering). MdnB's leader-recognition equivalent must be found by **structural superposition** of 5IG8 onto 5IG9, not by matching numbers. MdnB also binds the leader ~10× more weakly, so its recognition residues are expected to be less conserved.

### 5c. The two-tier framing *(the plan the focused run tested)*

- The **broad run** is a valid **superfamily-wide map** (universal core + DxR). Keep it.
- A **focused microviridin run** surfaces the **family-specific** layer (leader/substrate recognition).
- Comparing the two shows which residues are universal vs microviridin-specific. **This is now realized for MdnC (§5f).**

### 5d. The MdnC focused pipeline & the taxonomic-contamination discovery *(NEW — the key methodological result)*

The curation pipeline of §4c was executed on the RODEO **Group-1** ligase pool (**824** unique ester/amide ligase accessions).

**HMM classification** (`hmmscan --cut_ga`): **303 MdnC-type**, **215 MdnB-type**, **306 neither** (0 ambiguous). The bit-score delta between the two models is cleanly **bimodal** — the separation that pairwise identity never produced — confirming the HMM approach is the right classifier for this family.

**Two orthogonal validations of the classification:**
1. **RODEO gene-position (ligase-1 vs ligase-2 column) vs HMM call:** L1 is ~**87% MdnC-type**, L2 ~**85% MdnB-type** — a sequence-based method (HMM) independently reproducing a genomic-context signal (gene order), with **zero directional contradictions** (no sequence called ester by one and amide by the other).
2. **Old pairwise-identity classifier vs HMM:** again **zero directional contradictions**; the two methods disagree only on *stringency* (identity over-commits on borderline sequences the HMM abstains on), never on *direction* — mutual validation of the ester/amide split.

**The contamination discovery.** Inspecting the classified pool by phylum (from the RODEO spreadsheet) exposed the problem the broad run had only hinted at:
- The full 824-sequence Group-1 pool is **~58% non-cyanobacterial** (Bacteroidetes-dominated).
- The **303 MdnC-type set** breaks down as **172 Cyanobacteria + 117 Bacteroidetes + 12 Proteobacteria + 1 Acidobacteria + 1 uncultured** — i.e. **~43% non-cyanobacterial**.
- These non-cyanobacterial sequences are **genuine ester-forming graspetide ligases** (e.g. *Chryseobacterium*, *Sorangium*), so they legitimately clear the MvdD HMM — but they are **not microviridin**, and a small α7 alignment showed them **diverging at exactly residues 191/192/195** (Asp192→Lys, Asn195→Glu). Including them would recreate the broad-run contamination one taxonomic level finer.
- Crucially, the **475-bit threshold does not separate them by taxonomy** — it separates ester from amide chemistry; non-cyanobacterial ester-formers sit both above and below it (a coherent *Chryseobacterium* clade was found straddling the line). **Only a taxonomic filter removes them.**
- **Independent corroboration from BLAST:** a *wide* MdnC BLAST (down to ~43% identity) returned **~29% non-cyanobacterial** hits below ~68% identity, while a *focused* BLAST (≥69% identity) stayed **cleanly cyanobacterial**. Two different collection methods (RODEO enumeration and BLAST) independently show that identity cannot separate microviridin from distant ester-formers below a high identity floor.

**Set assembly after the taxonomic filter.** 172 cyanobacterial RODEO MdnC + **133** clean cyanobacterial focused-BLAST accessions (all 133 passed the MvdD HMM gate) → merged (279) → **CD-HIT 0.90 → 70 clusters**, + re-inserted anchor = **71 sequences**. (RODEO-alone would give 61 clusters; the BLAST set contributed **9 genuinely non-redundant, partly post-2021 sequences**, kept for the added evolutionary depth.)

**Phylogeny (IQ-TREE 2, model LG+I+R5 by BIC; logL −18163.8):**
- **Chamaesiphon minutus (WP_015158320)** — a deep-branching cyanobacterial MdnC — was placed near the outgroup by the fast FastTree preview but **resolved cleanly into the cyanobacterial in-group by IQ-TREE**; the FastTree placement was **long-branch attraction**. Its α7 signature (invariant Asp192, canonical Asn195, Ala191) independently confirms it as a genuine microviridin MdnC. Kept.
- Anchor **WP_012265553** sits with **AKE65548** (*Microcystis*) at **99/100** support, central in the clade.
- **Outgroup behaviour (a clean methodological footnote):** the 6 Bacteroidetes outgroup formed a proper monophyletic clan and rooted well, but the 2 Proteobacteria outgroup (CAN97345.1 / WP_012239784.1, identical) were long-branch-attracted into the cyanobacterial clade with weak support (UFBoot 70) and **both failed IQ-TREE's composition χ² test** (compositionally biased). This is an **outgroup artifact only** and does not affect the ConSurf set (all outgroup sequences are excluded by identity). *Lesson for future tree work: vet outgroup sequences for composition bias.*

### 5e. The focused MdnC ConSurf result *(NEW — the primary deliverable)*

Run on the **71-sequence** curated cyanobacterial MdnC MSA (external MSA; ConSurf's own homolog search bypassed), query WP_012265553 → 5IG9, Bayesian + LG. **All three target columns had 100% occupancy (71/71) and none was flagged "insufficient data"** — the grades are reliable.

| Residue | Grade | Confidence interval | Residue variety | Reading |
|---|---|---|---|---|
| **Asp192** | **9** (max), ConSurf-labeled **`f` = functional** | 9,9 (tight) | **D 100%** | Invariant; the standout result |
| **Asn195** | **5** (average) | 6,5 | N 71%, H 22% | Conserved but genuinely His-tolerant (a secondary contact) |
| **Glu191** | **4** (variable-leaning) | 4,4 | E 56%, D 26%, A 9% | Identity-variable but **~82% acidic** — charge, not identity, is conserved |

**Why the graded pattern is a stronger result than uniform conservation:**
- **Asp192 is a lone conservation spike in a variable loop** (neighbors: K190 = grade 2, E191 = 4, L193 = 8, D194 = 3, N195 = 5). A single invariant residue against a variable background is the textbook signature of a **specificity-determining position** held by specific functional selection — much more informative than a uniformly conserved region, where "essential" cannot be distinguished from "structurally rigid."
- **Glu191's low grade reflects identity variation, not lack of importance.** ConSurf scores identity; E191 shuffles among E/D/A while keeping the negative charge ~82% of the time. The functional constraint (the electrostatic contact with MdnA Arg17) is on the **charge**, which an identity-based method under-reports. A physicochemical-conservation view would show more constraint.
- **Asn195 is genuinely the most tolerant** (a real secondary contact).

**Independent convergence with the biochemistry.** The ConSurf ordering **D192 ≫ E191 ≈ N195** recapitulates Li 2016's mutagenesis: E191**K**/D192**K** (charge reversal) *abolished* activity, while E191**A**/D192**A** (milder) *retained partial activity*. So biochemistry says D192 is essential and E191 tolerates milder change — and evolution says exactly the same. Two fully independent lines (deep evolutionary selection and benchtop mutagenesis) converge on the same ranking. *Chamaesiphon*'s natural **Ala191 + Asp192** is a living instance of the "one acidic contact softened, the other intact" state the mutagenesis predicted to remain functional.

**Run-wide sanity:** 106 of 324 positions are grade 9 — a strong conserved core, consistent with a real ATP-grasp active site.

### 5f. Broad vs focused — the two-tier comparison, realized

| α7 residue | Broad run (§5a) | Focused run (§5e) | Change |
|---|---|---|---|
| Glu191 | grade 3 | grade 4 | ~flat (genuinely variable even within microviridin) |
| **Asp192** | grade 4 | **grade 9** | **+5 — the microviridin-specific signal surfaced** |
| Asn195 | grade 5 | grade 5 | flat (tolerant in both) |

This is the two-tier thesis confirmed on the actual target residue: **the focused, decontaminated set lifted Asp192 from variable to maximally conserved**, exactly the family-specific layer the broad superfamily map could not reveal. The residues that stayed low (E191, N195) did so because they are truly tolerant within microviridin, not because of set contamination — a distinction only the clean set could establish.

---

## 6. Decisions log & current state

**Decisions made (and why):**
- **Analyze MdnB and MdnC separately** — divergent clades; pooling averages away distinguishing residues.
- **Do not use the lab's private anti-SMASH deep-sea database** — unsorted, function-undetermined, biased, large; public databases are cleaner and defensible.
- **Precursor (MdnA) is not analyzed with ConSurf** — short, hypervariable, multicore; handle later via targeted alignment and/or ESM.
- **Fix the too-broad set by curating a microviridin-focused set**, not by nudging ConSurf's identity slider.
- **HMM family classification (TIGR04184/TIGR04185, `--cut_ga`) is the correct classifier**, not pairwise identity — the identity distribution is non-bimodal; the HMM split is clean and doubly validated.
- **A taxonomic (cyanobacteria-only) filter is required, not optional** — the decisive discovery of §5d. HMM chemistry classification alone leaves ~43% non-cyanobacterial contamination that diverges at the target residues.
- **Define the final set as a phylogenetic clade, confirmed by tree** — with the caveat that when a sequence is functionally validated by other means (e.g. *Chamaesiphon*'s α7 signature), that evidence overrides a shaky tree placement.
- **Upload the curated set to ConSurf as an external MSA** — to bypass ConSurf's own homolog search, which would re-contaminate.
- **Keep the BLAST-window set in the merge** — it added 9 non-redundant, partly post-2021 sequences; provenance trade-off noted (BLAST sequences are HMM+taxonomy-verified rather than RODEO-group-labeled).

**Current state:**
- **MdnC focused pipeline: COMPLETE through ConSurf.** Clean 71-sequence set built and validated; focused ConSurf run done; primary result in hand (§5e). Colored outputs (`consurf_colored_seq.pdf`, `5IG9_With_Conservation_Scores.pdb`) available for figures.
- The result is documented, reliable (100% target-column occupancy, no insufficient-data flags), and mechanistically coherent with Li 2016.

---

## 7. Homolog curation approach (detail) — *as actually executed*

This supersedes the original "BLAST-window first pass" plan. The executed pipeline (see §4c for rationale, §5d for results):

1. **RODEO backbone.** Take the Ramesh et al. 2021 graspetide dataset, pull **Group 1** ligase accessions (824 unique across the ligase-1/ligase-2 columns), carrying the spreadsheet's taxonomy columns.
2. **HMM classification.** `hmmscan --cut_ga` against combined TIGR04184 (`ATPgraspMvdD`) + TIGR04185 (`ATPgraspMvdC`); keep **MvdD-only** = MdnC-type (303).
3. **Taxonomic filter.** Keep `Phylum == Cyanobacteria` → **172** (this removes the 131 non-cyanobacterial ester-formers — the essential decontamination step).
4. **BLAST-window supplement.** Focused MdnC BLAST vs nr, 35–90% window, ≥70% coverage, 280–360 aa → **133** cyanobacterial accessions, each re-verified through the MvdD HMM gate.
5. **Merge + dedup.** Union (279) → **CD-HIT 0.90** → 70 clusters → **re-insert anchor** → **71**.
6. **Align + tree.** MAFFT L-INS-i (71 + 8 proximal non-cyano outgroup); trimAl `-gt 0.2` for the tree only; FastTree preview → IQ-TREE 2 (LG+I+R5). Confirm cyanobacterial in-group monophyly and anchor placement; note and set aside outgroup artifacts.
7. **ConSurf MSA.** Untrimmed alignment, clade only (outgroup removed), anchor included → readiness diagnostic (size 71, anchor present, 100% occupancy at 191/192/195, no leaks) → **submit as external MSA, LG model, Bayesian**.

**Key principle established:** *identity thresholds alone cannot define functional family membership.* HMM classification + taxonomic filtering + phylogenetic clade definition are all required; each addresses a distinct failure mode (chemistry, lineage, and residual outliers respectively).

---

## 8. Files & artifacts inventory

- **Literature review** — full background with verified citations. *(standalone file)*
- **ESM-2 methodology briefing** — PLM methods and plan. *(standalone file)*
- **ConSurf homolog-curation methodology review** — best practices for set composition, the taxonomic/clade rationale, and the graded-α7 grounding against Li 2016. *(standalone file — `MdnC_ConSurf_methodology_review_corrected.md`)*
- **This master document** — consolidated project reference.
- **Notebooks:** `MdnC_ConSurf_curation_pipeline.ipynb` — the consolidated pipeline (HMM classification B–C, taxonomic filter D, merge+CD-HIT E, align/trim/tree/readiness F). Week-2 and BLAST-window notebooks superseded.
- **HMM database:** `mvd_families.hmm` (+ `.h3m/.h3i/.h3f/.h3p`) — combined TIGR04184 + TIGR04185.
- **Curated sets:** `cyano_mdnc.fasta` (172), `blast_window.fasta` (133), `microviridin_mdnc_nr90.fasta` (71), `consurf_msa.fasta` (71, ConSurf input).
- **BLAST artifacts (raw XMLs removable):** `mdnC_query.fasta`; `blast_focused_hits.csv` + `blast_wide_hits.csv`; `blast_focused_window_cyano_accessions.txt` (133 — rebuild the BLAST set by efetch, no re-BLAST); `blast_wide_window_noncyano_accessions.txt` (250 — outgroup candidates).
- **Tree:** `iqtree_run.*` (treefile, iqtree report, contree, log). Model LG+I+R5.
- **ConSurf outputs (focused MdnC):** `5IG9_consurf_grades.txt`, `5IG9_With_Conservation_Scores.pdb`, `consurf_colored_seq.pdf`, `TheTree.txt`, `msa_aa_variety_percentage.csv`. **Plus** the earlier MdnC and MdnB *broad*-run grades + MSAs (keep as the superfamily-wide map).
- **Structures:** PDB 5IG9 (MdnC + MdnA leader), 5IG8 (MdnB), local copies in `data/`.
- **Lab-provided sequences:** fused MdnBC cyclase; multicore MdnA precursor (§2).
- **GitHub repo** — the working project folder, synced via GitHub Desktop.

---

## 9. Open questions & next steps

1. ~~Does the focused MdnC set fix α7?~~ **Done (§5e).** Asp192 rose to grade 9; E191/N195 graded lower for informative reasons.
2. **Two-run sensitivity analysis (optional, publication-grade):** re-run focused ConSurf with the divergent *Chryseobacterium* clade added and quantify how much the α7 grades drop — turns "we excluded contamination" into "we demonstrated microviridin-specificity."
3. **Physicochemical-conservation view of E191** — report an acidic-character (D/E-grouped) conservation alongside the raw grade, since identity-based grading under-reports its constraint.
4. **MdnB focused pipeline** — repeat §4c/§7 for MdnB (MvdC-only, cyanobacterial). Identify MdnB's leader-recognition residues by structural superposition of 5IG8 onto 5IG9 (not by number). Expect weaker conservation (MdnB binds the leader ~10× weaker).
5. **ESM-2 implementation** — per-residue and per-variant tracks; validate against §3; test whether ESM ranks D192 ≫ E191/N195 (a third independent line).
6. **Precursor (MdnA)** — separate sub-analysis; handle the multicore repeats carefully.
7. **Fused MdnBC handling** — score full-length vs per-domain; the H573E/H193E mutants are the eventual mechanistic ground truth.
8. **Colored-structure figures** — render `5IG9_With_Conservation_Scores.pdb` in ChimeraX with D192 highlighted at the MdnA-leader interface.

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

**Tools used in the focused pipeline (add full citations as needed):** HMMER (`hmmscan`/`hmmpress`); TIGRFAM/NCBIfam profiles TIGR04184, TIGR04185; MAFFT L-INS-i; trimAl; FastTree; IQ-TREE 2 with ModelFinder and UFBoot; BLASTP 2.17.0+.

*(Full annotated citations and a glossary are in the standalone literature-review, ESM-methodology, and ConSurf-methodology-review-corrected files.)*
