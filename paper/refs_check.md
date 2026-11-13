# refs.bib verification log

Built 2026-10-27. Starting point: the research notes, the lit-review report and `campaign/theory/LITERATURE.md`. Every entry was then checked against an external record:

- the Crossref REST API, one record per DOI;
- the arXiv API, for 20 ids;
- zbMATH Open, OpenAlex, the SMF page for Astérisque, and the Math-Net.ru table of contents, used where Crossref is silent.

No e-mail address or personal identifier was sent in any request.

## Summary

| Status | Count | Keys |
|---|---|---|
| Crossref DOI match (all metadata fields agree with the record) | 70 | every requested key except the three below |
| of which also cross-checked against an arXiv record | 18 | azizi2020, encinas2018, rodrigues2015, lee2011, janssen2005, lambiotte2008, russell2011, cox2013, cox2025, chatterjee2013, grassberger2006, villamartin2014, vazquez2010, llabres2026, eppstein2010, droz2003b, lavrentovich2014, park2015 |
| No DOI exists; verified in zbMATH (+ Math-Net.ru) | 2 | toom1980, kolmogorov1937 |
| arXiv only (confirmed unpublished: no journal_ref, no Crossref match) | 1 | hinrichsen2000fo |
| **UNVERIFIED (notes only)** | **0** | none; no `note = {UNVERIFIED: ...}` fields were needed |
| Extra entry, not in the requested list | 1 | lambiotte2007 (see the first item below) |

A field-by-field script compared refs.bib against the saved Crossref JSON: DOI, volume, first page or article number, year, author surnames and title. The only differences are the five explained below: henkel2008's subtitle, Crossref's mangled "Schl?gl", Crossref's splitting of "von der Twer" and "Arial Sznaider", and the Astérisque registration date.

## Needs your decision

1. **lambiotte2008: the arXiv id in the request points to a different paper.**
   - arXiv 0710.0914 is Lambiotte & Redner, "Dynamics of vacillating voters", J. Stat. Mech. (2007) L10001. It is a 2D square-lattice study, and it is the paper the lit-review report cites for a noise-free 2D active phase ("Lambiotte & Redner 2007").
   - The EPL paper "Dynamics of non-conservative voters" (EPL 82, 18007, 2008) is arXiv 0712.0364. The notes classify it as 1D only.
   - `lambiotte2008` keeps the EPL data, because the title, journal, volume, page, year and key all agree on it.
   - A verified extra entry, `lambiotte2007`, holds the vacillating-voters paper. If the sentence in the manuscript is about a 2D active phase, cite `lambiotte2007` instead.
2. **lipowsky1984: I chose Lipowsky, Kroll & Zia, PRB 27, 4499 (1983).**
   - Carlon et al. 1999 (Ref. 4, checked in the arXiv full text) cite exactly this paper for the "general, phenomenological" ω = 1/3 law, derived from an effective interface Hamiltonian. That is the same context as carlon1999 and selke1983.
   - It was published in **1983**. The key says 1984; I kept the key unchanged.
   - Companion option: R. Lipowsky, "Critical effects at complete wetting", Phys. Rev. B 32, 1731–1750 (1985), doi:10.1103/PhysRevB.32.1731 (verified in Crossref). It is the dedicated complete-wetting paper, with mean-field and fluctuation regimes. Add it if you want a complete-wetting source as such.
   - Fisher 1984 (fisher1984) is the random-walk/2D counterpart.

## Corrections relative to the requested descriptions

- **wang2015**: published in **J. Chem. Phys. 142, 164105 (2015)**, not PRE. Title: "Discontinuous non-equilibrium phase transition in a threshold Schloegl model for autocatalysis: Generic two-phase coexistence and metastability". Authors: C.-J. Wang, D.-J. Liu and J. W. Evans.
- **lavrentovich2014**: the first author is **M. O.** Lavrentovich (Maxim O.), not "O. Lavrentovich". Published as PRL 112, 138102 (2014), "Asymmetric mutualism in two- and three-dimensional range expansions". arXiv 1309.0273 matches.
- **carlon1999**: the title is "Interfacial adsorption in **two-dimensional** Potts models". "Planar Potts models" is the title of selke1983.
- **hinrichsen1997wet**: the title is "Model for nonequilibrium wetting transitions in two dimensions", with no leading "A".
- **eppstein2010**: pages are **71–97** per the Springer/Crossref chapter record. The notes and arXiv journal_ref say 71–98; I used the publisher's range.
- **encinas2018**: published as Sci. Rep. 8, 9338 (2018), titled "Fundamental ingredients for discontinuous phase transitions in the inertial majority vote model". The arXiv title differs ("...for the emergence of discontinuous phase transitions in the majority vote model"); the arXiv record's DOI field points to the Sci. Rep. paper.
- **azizi2020**: the published title reads "in **the** presence of" (the arXiv title omits "the").

## Items resolved from partial descriptions

- **droz2003b**: M. Droz and A. Lipowski, "Some recent developments in models with absorbing states", Braz. J. Phys. 33, 526–532 (2003). arXiv cond-mat/0308419 (journal_ref and DOI match).
- **rodrigues2015**: J. Stat. Mech. (2015) P01035, "Critical behavior in lattice models with two symmetric absorbing states". The guess was right.
- **lee2011**: M. Y. Lee and T. Vojta, "Generalized contact process with two symmetric absorbing states in two dimensions", PRE 83, 011114 (2011).
- **janssen2005**: H.-K. Janssen, "Survival and percolation probabilities in the field theory of growth models", J. Phys.: Condens. Matter 17, S1973–S1993 (2005). The arXiv record has no journal_ref; I matched it by title and author.
- **chatterjee2013**: S. Chatterjee and R. Durrett, "A first order phase transition in the threshold θ ≥ 2 contact process on random r-regular graphs and r-trees", Stoch. Process. Appl. 123, 561–578 (2013). The arXiv journal_ref agrees.
- **park2015**: S.-C. Park (single author), "Order–disorder transition in the two-dimensional interacting monomer-dimer model: Ising criticality", J. Stat. Mech. (2015) P10009. arXiv 1507.05537.
- **monetti1995**: confirmed as PRE 52, 5825–5831 (1995), "Critical edge between frozen extinction and chaotic life".
- **villamartin2015**: confirmed as PNAS 112(15), E1828–E1836 (2015). Crossref has no page, so the pages come from OpenAlex.
- **menyhard1994** and **bramson1980**: Crossref bibliographic search found the DOIs 10.1088/0305-4470/27/18/022 and 10.1007/BF01013315.
- **cox2013**: Astérisque 349 (2013), vi+113 pp., doi:10.24033/ast.921. Crossref's 2018 date is the DOI registration date; the SMF page gives the volume and year. No page range is given because the paper fills the whole volume.

## Style choices (so you can override them)

- apsrev4-2 sentence-cases article titles: `format.title` applies `"t" change.case$`, which I checked in the REVTeX 2023 bst. Proper nouns and math are therefore braced, for example {I}sing, {P}otts, {L}otka--{V}olterra, {S}chl{\"o}gl/{S}chloegl, {L}angevin, {G}ibbs, {L}ife, {$q$} and {$\mathbb{Z}^d$}. main.tex loads amssymb, so `\mathbb` is safe.
- "Voter" is braced as {V}oter only where the published title capitalizes it (dornic2001, holley1975, cox1986), following your example list. Elsewhere "voter" stays lowercase, as published. Remove those braces if you prefer "voter model" throughout.
- Avrami's title is protected verbatim ("Kinetics of Phase Change. I General Theory"), because sentence-casing would mangle the "I".
- Journal names use APS/ISO4 abbreviations. Article numbers go in `pages`.
- eprint/archivePrefix appears only on the arXiv-only item, because apsrev4-2 would also print eprints for published papers.
- Author names are given as initials. The bylines vary, so I normalized them:
  - Cox and Perkins appear as J. T. Cox and E. A. Perkins in all four Cox entries (cox2025's byline says "Ted Cox, Edwin Perkins").
  - Janssen is "H. K." in 1981 and "H.-K." in 2005, as in each byline; it is the same person.
  - Sznaider is "G. A."; see szabo2004 below.
- I could not compile, because there is no TeX installation on this machine. Instead, a Python parser checked brace balance, required fields per entry type, ASCII-only content and key coverage: all 73 requested keys are present, plus lambiotte2007, and all 47 keys currently cited in main.tex resolve.

## Per-key log

| Key | Verified via | Notes / flags |
|---|---|---|
| hinrichsen2000 | Crossref DOI match | Adv. Phys. 49(7), 815–958 (2000) |
| odor2004 | Crossref DOI match | RMP 76, 663–724 (2004) |
| henkel2008 | Crossref DOI match (book) + zbMATH 1165.82002 + OpenAlex | Crossref lacks the authors and the "Vol. 1: Absorbing Phase Transitions" subtitle; both come from zbMATH/OpenAlex. Springer, Dordrecht |
| lubeck2004 | Crossref DOI match | IJMPB 18(31–32), 3977–4118 (2004) |
| janssen1981 | Crossref DOI match | Title: "On the nonequilibrium phase transition in reaction-diffusion systems with an absorbing stationary state" |
| grassberger1982 | Crossref DOI match | "On phase transitions in Schlögl's second model" (Crossref garbles the ö) |
| grassberger1984 | Crossref DOI match | "A new type of kinetic critical phenomenon"; L105–L109 |
| menyhard1994 | Crossref (DOI found by search) | "One-dimensional non-equilibrium kinetic Ising models with branching annihilating random walk"; Crossref omits the accent in Menyhárd |
| hinrichsen1997 | Crossref DOI match | PRE 55, 219–226 (1997) |
| dornic2001 | Crossref DOI match | PRL 87, 045701 (2001) |
| alhammal2005 | Crossref DOI match | PRL 94, 230601 (2005) |
| vazquez2008 | Crossref DOI match | PRE 78, 061127 (2008) |
| droz2003 | Crossref DOI match | PRE 67, 056108 (2003) |
| lipowski2002 | Crossref DOI match | PRE 65, 056114 (2002) |
| droz2003b | Crossref DOI match + arXiv cond-mat/0308419 | Resolved to Braz. J. Phys. 33, 526 (2003); see above |
| castellano2009 | Crossref DOI match | PRE 80, 041129 (2009) |
| azizi2020 | Crossref DOI match + arXiv 2007.10933 | Published title has "in the presence of" |
| rodrigues2015 | Crossref DOI match + arXiv 1501.01116 | J. Stat. Mech. (2015) P01035 |
| lee2011 | Crossref DOI match + arXiv 1010.3298 | PRE 83, 011114 (2011) |
| canet2005 | Crossref DOI match | "Nonperturbative fixed point in a nonequilibrium phase transition" |
| janssen2005 | Crossref DOI match + arXiv cond-mat/0304631 | J. Phys.: Condens. Matter 17, S1973 (2005) |
| molofsky1999 | Crossref DOI match | Theor. Popul. Biol. 55(3), 270–282 (1999) |
| lambiotte2008 | Crossref DOI match + arXiv 0712.0364 | **The arXiv id in the request (0710.0914) is a different paper; see "Needs your decision"** |
| russell2011 | Crossref DOI match + arXiv 1012.1988 | PRL 106, 165702 (2011) |
| encinas2018 | Crossref DOI match + arXiv 1802.05256 | Sci. Rep. 8, 9338 (2018); published title differs from the arXiv title |
| clifford1973 | Crossref DOI match | "A model for spatial conflict", Biometrika 60(3), 581–588 |
| holley1975 | Crossref DOI match | "Ergodic theorems for weakly interacting infinite systems and the voter model", 643–663 |
| liggett1985 | Crossref DOI match + zbMATH 0559.60078 | Grundlehren der mathematischen Wissenschaften vol. 276 (number from zbMATH) |
| cox1986 | Crossref DOI match | Ann. Probab. 14(2), 347–370 (1986) |
| bramson1980 | Crossref (DOI found by search) | Z. Wahrsch. verw. Gebiete 53(2), 183–196 (1980) |
| krapivsky1992 | Crossref DOI match | PRA 45, 1067–1072 (1992) |
| cox2013 | Crossref DOI match + SMF page + arXiv 1103.1676 | Volume 349 and year 2013 from SMF; no pages (whole volume) |
| cox2008 | Crossref DOI match | Ann. Appl. Probab. 18(2), 747–812 (2008) |
| cox2010 | Crossref DOI match | EJP 15, 1190–1266 (2010) |
| cox2025 | Crossref DOI match + arXiv 2310.00812 | EJP 30, pp. 1–86 (2025). The EJP "paper no." was not checked (Project Euclid blocks automated access), so apsrev4-2 will print "30, 1 (2025)". Minor |
| chatterjee2013 | Crossref (DOI found by search) + arXiv 1011.1567 | SPA 123(2), 561–578 (2013) |
| bezuidenhout1994 | Crossref DOI match | Ann. Probab. 22(3), 1160–1194 (1994) |
| ziff1986 | Crossref DOI match | PRL 56, 2553–2556 (1986) |
| bidaux1989 | Crossref DOI match | "Order of the transition versus space dimension in a family of cellular automata" |
| lubeck2006 | Crossref DOI match | J. Stat. Phys. 123(1), 193–221 (2006) |
| grassberger2006 | Crossref DOI match + arXiv cond-mat/0510428 | J. Stat. Mech. (2006) P01004 |
| villamartin2014 | Crossref DOI match + arXiv 1402.1406 | PRE 89, 012145 (2014) |
| villamartin2015 | Crossref DOI match + OpenAlex | Pages E1828–E1836 come from OpenAlex (Crossref has none) |
| hinrichsen2000fo | arXiv cond-mat/0006212 | arXiv only (v1, 13 Jun 2000). No journal_ref, and a Crossref title search found no journal version |
| evans1991 | Crossref DOI match | "Characterizing kinetics near a first-order catalytic-poisoning transition" |
| szabo2004 | Crossref DOI match | Crossref splits the name as "Arial Sznaider, Gustavo"; Crossref has the same author as "Sznaider, Gustavo Ariel" (PRE 76, 051921), so the entry uses "Sznaider, G. A." |
| szabo2007 | Crossref DOI match | Phys. Rep. 446(4–6), 97–216 (2007) |
| guisoni2013 | Crossref DOI match | "Phase diagram of a cyclic predator-prey model with neutral-pair exchange" |
| vazquez2010 | Crossref DOI match + arXiv 1002.1251 | J. Stat. Mech. (2010) P04007 |
| llabres2026 | Crossref DOI match (has-preprint relation to arXiv 2505.11358) + arXiv | PRR 8, 013015 (2026), doi:10.1103/6xdk-2cjy |
| bennett1985 | Crossref DOI match | PRL 55, 657–660 (1985) |
| toom1980 | zbMATH Zbl 0441.68053 (chapter) + Zbl 0425.00018 (volume) | No DOI exists. In *Multicomponent Random Systems*, ed. by R. L. Dobrushin and Ya. G. Sinai (D. Griffeath contributing editor), Adv. Probab. Relat. Top. 6 (Dekker, New York, 1980), pp. 549–575 |
| wang2015 | Crossref DOI match | **Journal corrected: J. Chem. Phys. 142, 164105 (2015), not PRE** |
| liu2018 | Crossref DOI match | PRL 121, 120603 (2018) |
| liu2025 | Crossref DOI match | PRE 111, 064134 (2025), doi:10.1103/y54y-5dvx |
| rikvold1994 | Crossref DOI match | PRE 49, 5080–5090 (1994) |
| avrami1939 | Crossref DOI match | J. Chem. Phys. 7, 1103–1112 (1939); title protected verbatim |
| kolmogorov1937 | zbMATH/JFM 63.1130.03 + Math-Net.ru issue contents | No DOI. Izv. Akad. Nauk SSSR Ser. Mat. 1(3), 355–359 (1937). The title is the indexed German one, "Zur Statistik der Kristallisationsvorgänge in Metallen". It is often cited in English as "On the statistical theory of the crystallization of metals"; swap it in if you prefer |
| carlon1999 | Crossref DOI match | **Title is "...two-dimensional Potts models"**; J. Stat. Phys. 96, 531–543 (1999) |
| selke1983 | Crossref DOI match | Z. Phys. B 50, 113–116 (1983) |
| fisher1984 | Crossref DOI match | J. Stat. Phys. 34, 667–729 (1984) |
| lipowsky1984 | Crossref DOI match + Carlon et al. reference list | **PRB 27, 4499 (1983): year 1983, not 1984**; see "Needs your decision" |
| hinrichsen1997wet | Crossref DOI match | Title has no leading "A"; PRL 79, 2710–2713 (1997) |
| park2015 | Crossref DOI match + arXiv 1507.05537 | S.-C. Park, single author |
| deoliveira2015 | Crossref DOI match | PRE 92, 062126 (2015) |
| eppstein2010 | Crossref DOI match (chapter + book record) + arXiv 0911.2890 | Ed. A. Adamatzky (Springer, London). Pages 71–97 (publisher) vs 71–98 (arXiv/notes) |
| blok1999 | Crossref DOI match | PRE 59, 3876–3879 (1999) |
| fates2005 | Crossref DOI match (10.25088/ComplexSystems.16.1.1) | Crossref's author list has two junk entries (e-mail/affiliation strings). The real authors are Nazim A. Fatès and Michel Morvan. Complex Syst. 16(1), 1–27 |
| munoz1996 | Crossref DOI match | "Critical behavior of systems with many absorbing states", PRL 76, 451–454 |
| salmon2011 | Crossref DOI match + OpenAlex | SC '11 proceedings (ACM), pp. 1–12 |
| monetti1995 | Crossref DOI match | Confirmed PRE 52, 5825–5831 (1995) |
| lavrentovich2014 | Crossref DOI match + arXiv 1309.0273 | **First author M. O. Lavrentovich**; PRL 112, 138102 (2014) |
| kpz1986 | Crossref DOI match | PRL 56, 889–892 (1986) |
| lambiotte2007 (extra) | Crossref DOI match + arXiv 0710.0914 | "Dynamics of vacillating voters", J. Stat. Mech. (2007) L10001. Not cited yet |

## Additions 2026-10-30 (referee round 3: literature and novelty)

Twenty-five entries were added for the rewritten Introduction, Sec. IV D, Sec. VI and the Discussion. Each was looked up through the Crossref REST API (`query.bibliographic`, no personal identifier sent); the DOI, journal, volume, pages and year below come from the returned record. Guisoni et al. 2013 was already present; its full text has now been read (CONICET open-access copy).

| Key | Record | Notes |
|---|---|---|
| cox1991 | 10.1007/978-1-4612-0459-6_8, pp. 189–201 (1991) | Chapter in the Spitzer Festschrift; editors and series from the book record (Durrett & Kesten, Progress in Probability 28): check the printed volume |
| liggett1994 | 10.1214/aop/1176988729, Ann. Probab. 22, 764–802 | A correction exists: Ann. Probab. 27, 2120 (1999) |
| durrett1994 | 10.1098/rstb.1994.0028, Phil. Trans. R. Soc. B 343, 329–350 | |
| durrett1994b | 10.1214/aop/1176988861, Ann. Probab. 22, 289–333 | Durrett & Neuhauser |
| durrett2009 | 10.1214/08-AAP590, Ann. Appl. Probab. 19, 477–496 | |
| schweitzer2009 | 10.1140/epjb/e2009-00001-3, EPJ B 67, 301–318 | |
| bagnoli2001 | 10.1103/PhysRevE.63.046116 | |
| bagnoli2005 | 10.1103/PhysRevE.71.046108 | |
| bassler1996 | 10.1103/PhysRevLett.77.4094, PRL 77, 4094–4097 | Cited only through Bagnoli et al.; not read |
| azizi2018 | 10.1088/1742-5468/aae2dd, J. Stat. Mech. (2018) 103208 | |
| szabo2005 | 10.1088/0305-4470/38/30/005, J. Phys. A 38, 6689–6702 | Source of µ_cr = 0.02662(2) quoted in Szabó & Fáth |
| szabo2007b | 10.1103/PhysRevE.76.051921 | Szabó, Szolnoki, Sznaider; Ref. [16] of Guisoni et al. |
| szabo2008 | 10.1103/PhysRevE.77.011906 | Title completed from memory beyond the Crossref snippet ("...with four or six species"): verify |
| barghathi2016 | 10.1103/PhysRevE.93.022120 | |
| deoliveira1993 | 10.1088/0305-4470/26/10/006, J. Phys. A 26, 2317–2324 | |
| schulman1978 | 10.1007/BF01011727, J. Stat. Phys. 19, 293–314 | |
| monetti2001 | 10.1103/PhysRevE.65.016103 | |
| adachi2004 | 10.1016/j.physd.2004.04.010, Physica D 198, 182–196 | |
| fates2010 | 10.1007/978-1-84996-217-9_14, pp. 257–274 | |
| roy2024 | 10.1007/978-3-031-56943-2_15, CCIS 2021, pp. 192–203 | Abstract only; editors from the Springer page |
| odor1996 | 10.1103/PhysRevE.53.2231, PRE 53, 2231–2238 | |
| tome1993 | 10.1103/PhysRevE.47.948, PRE 47, 948–952 | |
| munoz1998 | 10.1209/epl/i1998-00123-2, EPL 41, 147–152 | |
| delossantos2003 | 10.1103/PhysRevE.67.021607 | |
| munoz2005 | 10.1140/epjb/e2005-00029-3, EPJ B 43, 73–79 | |

Other changes: `lipowsky1984` renamed `lipowsky1983` (the paper is from 1983; the citation in main.tex was updated). `lambiotte2008` (the 1D EPL paper, never cited) removed; `lambiotte2007` stays.
