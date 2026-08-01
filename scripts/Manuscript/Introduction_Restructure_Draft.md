# Introduction Restructure — Sentence Pull

Sentences pulled from the current Introduction (SDManuscript.Rmd, lines 55-108) and sorted into the six planned paragraphs. Citations preserved as-is. A few sentences could plausibly sit in more than one bucket — these are flagged with *[also fits: ...]*.

---

## P1 — General intro: biodiversity crisis, use of ENMs
*(first half currently pretty good — mostly draws from current paragraph 1)*

- The effects of anthropogenic stressors have led to a global extinction crisis, with estimates suggesting that 20-45% of plant species are facing extinction [@bachman2024extinct; @brummitt2015green; @pimm2015many; @nic2020extinction].
- Determining which plant species to focus our conservation efforts (e.g., active restoration, preserve creation, *ex situ* collections) on requires an array of natural history data seldom available to decision-makers [@heywood2017plant].
- The environmental niche of a species, the geographic extent of this niche that is occupied, the occupancy of sites within the range, and the census size of populations are all essential attributes for evaluating the species against anthropogenic stressors [@natural2001iucn; @faber2012natureserve; @usfws2016ssa].
- Although these attributes are relatively simple to characterise on paper, implementation is time-consuming, requiring extensive fieldwork and travel to field sites - therefore requiring proxies for use in conservation assessments [@juffe2016assessing; @bland2015cost; @pelletier2018predicting].
- Environmental niche models (ENMs or Species Distribution Models SDMs) have made enormous contributions for two of these attributes, describing the realised environmental niche of a species, and the species' geographic extent within it [@syfert2014using; @anderson2023integrating; @kass2021improving].
- Further, this progress - paired with a growing volume of occurrence data from citizen science platforms and the digitisation of natural history museum records - offers promise for modelling additional population parameters [@markham2023review; @feldman2021trends]. *[also fits: P4/P5 — sets up new-population and census-size modelling]*
- An ENM predicts a single variable: the probability of suitable habitat for the taxon (taddeo2026diversitydistributions). *[also fits: P4 — sets up the leap from suitability to realised distribution]*

---

## P2 — Field validation of models, across a species' range

- However, despite these advances, the utility of models from high-resolution data is seldom field verified, especially at landscape scales [@chiffard2020anbs; @a2022species].
- Collecting data on whether areas are favourable to continued recruitment is perhaps more pertinent than assessing whether long-lived individuals persist [taddeo2026diversitydistributions].
- ENMs generally suffer from having few, often spatially biased, occurrence records as dependent variables, which often fail to characterise the ecological breadth of the species [@stolar2015accounting; @feeley2011keep].
- Accordingly, while many ENMs achieve high performance on small hold-out subsets, they are unlikely to detect many new populations during ground verification [@a2022species].
- To increase the number of presences and absences in sites that are relatively similar to those harboring presences, which can be used for training models, adaptive-niche-based sampling (ANBS) has become increasingly employed, wherein cells with a high probability of occurrence or where models have higher prediction uncertainty are preferentially visited [@guisan2006using].
- Refitting a model with this newly acquired data not only allows for the acquisition of a larger number of presences and absences, but also allows for verifying coordinate placement, confirming that historic points are still extant, reducing model uncertainty, and acquiring additional data such as census estimates and life stages [@stockwell2002effects; @wisz2008effects; @jansen2022stop; @mondain2024adaptive]. *[also fits: P5 — mentions census estimates]*
- However, ground-truthing results are seldom shared or reported for plant species [@johnson2023field; @borokini2023iterative; @zimmer2023field].
- Accordingly, ground verification has not been conducted for them, nor has apparently much effort been made to improve ground verification for rare species... *[also fits: P4 — follows the dispersal/connectivity modelling discussion]*

---

## P3 — Effect of different-resolution predictors on model results

- However, the historic mismatch between the resolution of variables governing species distributions and the data available to serve as predictors of environmental conditions has restricted the interpretation and implementation of these models in highly heterogeneous environments [@guisan2013predicting].
- Recent advances in remote sensing, statistics, and software, have narrowed this resolution mismatch, providing environmental predictors which accurately reflect local environment even in highly heterogeneous environments [CITE].
- Accordingly, most of our knowledge about producing ENMs relies on simulated species and data, and generally at spatial resolutions useful for making large scale conservation policy decisions, while data at the resolution for a land manager to make individual decisions is nascent.
- An historic complication with the implementation and interpretation of ENMs is a mismatch between the spatial resolution of the independent variables available to model the species' fundamental niche and the factors governing the true distribution of populations - the realised niche [@carscadden2020niche; @chauvier2022resolution; @lembrechts2019incorporating].
- Recent papers have had mixed results regarding the effects of imprecisely mapped occurrence data on ENM predictions, with indications that models generated in more heterogeneous environments, and at finer resolutions (e.g. ca. 3 arc-second to 10 arc-minutes (~14.5km at 38$^\circ$)) suffer minor decreases in model performance [@graham2008influence; @smith2023including] with real species, while increasing error in mapping has drastic effects on model predictions with virtually simulated species [@gabor2022positional].
- A further mismatch of resolution is the year in which data on geographic localities were obtained and current conditions, which allow for positive population growth [@bracken2022maximizing].
- Historic occurrence data may now reflect conditions inhospitable to population maintenance, or even represent populations that were simply sinks from more robust populations [@bracken2022maximizing]; using these records may decrease model performance.
- High-resolution spatial data eventually enable the determination of the extent of individual populations by detecting the edges of suitable habitat. *[also fits: P5 — population boundary mapping]*

---

## P4 — Finding new populations and subpopulations using a model

- Rather in most conservation applications, the insight desired from an ENM is the species' realised distribution (cite).
- The bridge between an ENM and a plant population's presence is related to the dispersal of propagules and the establishment of the population, rather than the distribution of the fundamental niche alone [@engler2009migclim].
- Calls - and significant progress have been made for better integration of SDM predictions with models of dispersal since the rise in popularity of ENMs [@guisan2005predicting; @franklin2010moving].
- Most of this progress focuses on the goal of 'connecting' currently suitable habitats to habitats predicted to be suitable under climate change scenarios [@engler2009migclim; @bocedi2014range].
- To assist practitioners in detecting new populations or extending the range of known populations, we propose modelling plant occupancy as a random variable that depends on distance from known occupied sites.
- In simple terms, this perspective links habitat suitability to the tenets of island biogeography (CITE).
- 'Distance' may be defined as Euclidean (or Haversine for large distances), or as a least-cost distance reflecting a generalised surface which conveys the difficulty for seeds to travel between the nearest occupied sites and the site of interest [@etherington2016least].
- By combining habitat suitability with both Euclidean and least-cost distance, we distinguish candidate sites as extensions of known populations, new sub-populations, or new populations, rather than relying on suitability alone to prioritise field verification.

---

## P5 — Estimating population census size

- Obtaining reliable estimates of plant population census sizes is a time-intensive process requiring two pieces of information: 1) measurements of plant density and 2) population extent.
- Often, instead of directly determining either the population's boundaries or density, observers estimate population sizes; however, these estimates are often unreliable [@ReischChristoph2018Acom].
- Several promising field methods for acquiring density measurements have recently been developed [@rominger2019application; @krening2021sampling; @ermakova2021densities; @alfaro2019optimal; @schorr2013using] to estimate better effective population census size (*N*), as well as progress in the application of genomic methods to estimate effective population sizes (*N~e~*) (e.g. linkage disequilibrium) [@waples2024practical].
- However, these methods are focused on descriptive rather than predictive processes, estimating uncertainty in measurements *within* a population rather than across the species' range.
- We propose that estimates of plant density are best generated using methods that allow fitting statistical models incorporating spatial covariates and that predict estimates and measurements of uncertainty across gridded surfaces [@oliver2012population; @doser2022spabundance; @a2022species].
- Mapping the boundaries of a population is another time-consuming, albeit essential, task for realistic estimates of population census size.
- Currently, population boundaries are implicitly delineated by practitioners walking distances beyond which gene flow between individuals should be minimal (e.g., 1km), searching for additional individuals, or via aerial imagery [@alma9955306914202441].
- Many rare plant species are small and inconspicuous, and many taxa grow in small, clustered groups in specific habitats, making detection via either method difficult [@chen2013imperfect; @condit2000spatial].

---

## P6 — All hypotheses together
*(currently pretty good — this is the existing "Here we:" list, lines 104-107; carry forward largely as-is)*

- Here we: 1) evaluate the results of the best-performing model with groundtruthing (within 450m of existing occurrences).
- 2) Compare the results of models from different resolutions to the best model, as based on hold-out test data.
- 3) Determine whether model predictions are useful to find new sub-populations (450m - 990m from existing occurrences), or populations (> 990m from existing occurrences).
- 4) Estimate population census size using 4a) estimates of plant density, and 4b) modelling of population boundaries.

---

## Notes / gaps to fill during rewrite
- P2 (field validation) is currently thin on "across a species' range" — may need a new sentence about spatial coverage/representativeness of validation, not just ground-truthing volume.
- Several `[CITE]` / `(cite)` / `(taddeo2026diversitydistributions)` placeholders carried over — still need real citation keys before this becomes final text.
- The sentence "An ENM predicts a single variable..." and "Further, this progress... offers promise for modelling additional population parameters" both act as bridges — decide whether they open P1 (as a closer) or open P4/P5 as a topic sentence.
