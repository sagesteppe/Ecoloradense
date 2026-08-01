# Introduction Restructure — Sentence Pull

---

## P1 — General intro: biodiversity crisis, use of ENMs
*(first half currently pretty good — mostly draws from current paragraph 1)*

- The effects of anthropogenic stressors have led to a global extinction crisis, with estimates suggesting that 20-45% of plant species are facing extinction [@bachman2024extinct; @brummitt2015green; @pimm2015many; @nic2020extinction].
- Determining which plant species to focus our conservation efforts (e.g., active restoration, preserve creation, *ex situ* collections) on requires an array of natural history data seldom available to decision-makers [@heywood2017plant].
- The environmental niche of a species, the geographic extent of this niche that is occupied, the occupancy of sites within the range, and the census size of populations are all essential attributes for evaluating the species against anthropogenic stressors [@natural2001iucn; @faber2012natureserve; @usfws2016ssa].
- Although these attributes are relatively simple to characterise on paper, implementation is time-consuming, requiring extensive fieldwork and travel to field sites - therefore requiring proxies for use in conservation assessments [@juffe2016assessing; @bland2015cost; @pelletier2018predicting].
- Environmental niche models (ENMs or Species Distribution Models SDMs) have made enormous contributions for two of these attributes, describing the realised environmental niche of a species, and the species' geographic extent within it [@syfert2014using; @anderson2023integrating; @kass2021improving].
- Further, this progress - paired with a growing volume of occurrence data from citizen science platforms and the digitisation of natural history museum records - offers promise for modelling additional population parameters [@markham2023review; @feldman2021trends].

### Draft prose

The effects of anthropogenic stressors have led to a global extinction crisis, with estimates suggesting that 20-45% of plant species are facing extinction [@bachman2024extinct; @brummitt2015green; @pimm2015many; @nic2020extinction].
Determining which plant species to focus our conservation efforts (e.g., active restoration, preserve creation, *ex situ* collections) on requires an array of natural history data seldom available to decision-makers [@heywood2017plant].
The environmental niche of a species, its occupied geographic extent, the occupancy of sites within the range, and the census size of populations are all essential attributes for evaluating the species against anthropogenic stressors [@natural2001iucn; @faber2012natureserve; @usfws2016ssa].
Although these attributes are relatively simple to characterise on paper, implementation is time-consuming, requiring extensive fieldwork - therefore requiring proxies for use in conservation assessments [@juffe2016assessing; @bland2015cost; @pelletier2018predicting].
Environmental niche models (ENMs or Species Distribution Models, SDMs) have made enormous contributions to two of these attributes, describing the realised environmental niche of a species, and the species' geographic extent within it [@syfert2014using; @anderson2023integrating; @kass2021improving].
Further, this progress - paired with a growing volume of occurrence data from citizen science platforms and the digitisation of natural history museum records - offers promise for modelling the remaining two attributes, occupancy and census size [@markham2023review; @feldman2021trends].

*Changes from the source text: tightened "for two of these attributes" → "to two of these attributes"; reworked the closing clause from "offers promise for modelling additional population parameters" into an explicit callback to the four attributes named earlier (niche, extent, occupancy, census size) and a forward pointer into P2-P5, since that sentence now has to do bridging work it didn't before. Wordiness pass: "the geographic extent of this niche that is occupied" → "its occupied geographic extent"; dropped redundant "and travel to field sites" (fieldwork implies it).*

---

## P2 — Field validation of models, across a species' range

- However, despite this progress, the utility of models is seldom field verified, especially at landscape scales [@chiffard2020anbs; @a2022species].
- ENMs generally suffer from having few, often spatially biased, occurrence records as dependent variables, which often fail to characterise the ecological breadth of the species [@stolar2015accounting; @feeley2011keep].
- Accordingly, while many ENMs achieve high performance on small hold-out subsets, there might be biases in the training data carried onto field predictions [@a2022species].
- To increase the number of presences and absences in sites that are relatively similar to those harboring presences, which can be used for training models, adaptive-niche-based sampling (ANBS) has become increasingly employed, wherein cells with a high probability of occurrence or where models have higher prediction uncertainty are preferentially visited [@guisan2006using].
- Refitting a model with this newly acquired data not only allows for the acquisition of a larger number of presences and absences, but also allows for verifying coordinate placement, confirming that historic points are still extant, reducing model uncertainty, and acquiring additional data such as census estimates and life stage demographics [@stockwell2002effects; @wisz2008effects; @jansen2022stop; @mondain2024adaptive]. 
- However, ground-truthing results are seldom shared or reported for plant species [@johnson2023field; @borokini2023iterative; @zimmer2023field].
- Accordingly, ground verification has not been conducted for them, nor has apparently much effort been made to improve ground verification for rare species... 

### Draft prose

However, despite this progress, the utility of models is seldom field verified, especially at landscape scales [@chiffard2020anbs; @a2022species].
ENMs generally suffer from having few, often spatially biased, occurrence records as dependent variables, which often fail to characterise the ecological breadth of the species [@stolar2015accounting; @feeley2011keep].
Accordingly, while many ENMs achieve high performance on small hold-out subsets, there might be biases in the training data carried onto field predictions [@a2022species].
To generate more presence and absence records for training, adaptive-niche-based sampling (ANBS) has become increasingly employed, preferentially visiting cells with high predicted occurrence probability or high prediction uncertainty [@guisan2006using].
Refitting a model with this newly acquired data not only expands the presence/absence pool, but also verifies coordinate placement, confirms that historic points are still extant, reduces model uncertainty, and yields additional data such as census estimates and life stage demographics [@stockwell2002effects; @wisz2008effects; @jansen2022stop; @mondain2024adaptive].
However, ground-truthing results are seldom shared or reported for plant species, nor has much effort apparently been made to improve ground verification for rare species specifically [@johnson2023field; @borokini2023iterative; @zimmer2023field].

*Wordiness pass: tightened the ANBS sentence and the "not only allows for... but also allows for..." refitting sentence.*

**Still open (per the notes below):** this paragraph is heavy on field-validation-is-rare/ANBS-as-a-fix but doesn't yet say much about validation "across a species' range" specifically — may need a sentence on spatial coverage/representativeness, not just volume of ground-truthing.

---

## P3 — Effect of different-resolution predictors on model results
The spatial-grain of environmental variables have often been too coarse to model micro-environments which govern many rare plants distributions, especially in heterogenous environments [@guisan2013predicting].
- Accordingly, most of early ENM modelling has focused models useful for conservation policy decisions, while models at the resolution for a natural resource professionals to make individual management decisions is nascent.
- Recent advances in remote sensing, statistics, and software, have narrowed this resolution mismatch, [CITE].
- Recent papers have had mixed results regarding the effects of imprecisely mapped occurrence data on ENM predictions, with indications that models generated in more heterogeneous environments, and at finer resolutions (e.g. ca. 3 arc-second to 10 arc-minutes (~14.5km at 38$^\circ$)) suffer minor decreases in model performance [@graham2008influence; @smith2023including] with real species, while increasing error in mapping has drastic effects on model predictions with virtually simulated species [@gabor2022positional].

### Draft prose

The spatial grain of environmental variables has often been too coarse to model the micro-environments that govern many rare plants' distributions, especially in heterogeneous environments [@guisan2013predicting].
Accordingly, most early ENM work targeted resolutions useful for conservation policy decisions, while resolutions fine enough for natural resource professionals' individual management decisions remain nascent.
Recent advances in remote sensing, statistics, and software have narrowed this resolution mismatch [CITE].
Recent papers have had mixed results regarding the effects of imprecisely mapped occurrence data on ENM predictions, with indications that models generated in more heterogeneous environments, and at finer resolutions (e.g. ca. 3 arc-second to 10 arc-minutes (~14.5km at 38$^\circ$)) suffer minor decreases in model performance [@graham2008influence; @smith2023including] with real species, while increasing error in mapping has drastic effects on model predictions with virtually simulated species [@gabor2022positional].

*Wordiness pass: "focused on models useful for... models at the resolution needed for..." repeated "models" awkwardly — reworded to "targeted resolutions useful for... resolutions fine enough for...".*

---


## P4 — Finding new populations and subpopulations using a model
- An ENM predicts a single variable: the probability of suitable habitat for the taxon (taddeo2026diversitydistributions). 
- However, in most conservation applications, the insight desired from an ENM is the species' realised distribution (cite).
- Significant progress has been made integrating SDM predictions with models of dispersal, recognising that the bridge between an ENM and a plant population's presence is propagule dispersal and establishment, not the distribution of the fundamental niche alone [@guisan2005predicting; @franklin2010moving; @engler2009migclim].
- Previous innovations have largely focused on connecting currently suitable habitat to habitat predicted to be suitable under climate change [@engler2009migclim; @bocedi2014range].
- To assist practitioners in detecting new populations or extending the range of known populations, we propose modelling plant occupancy as a random variable that depends on distance from known occupied sites.
- 'Distance' may be defined as Euclidean (or Haversine for large distances), or as a least-cost distance reflecting a generalised surface which conveys the difficulty for seeds to travel between the nearest occupied sites and the site of interest [@etherington2016least].
- By combining habitat suitability with both Euclidean and least-cost distance, we distinguish candidate sites as extensions of known populations, new sub-populations, or new populations, rather than relying on suitability alone to prioritise field verification.

### Draft prose

An ENM predicts a single variable: the probability of suitable habitat for the taxon (taddeo2026diversitydistributions).
However, in most conservation applications, the insight desired from an ENM is the species' realised distribution (cite).
Significant progress has been made integrating SDM predictions with models of dispersal, recognising that the bridge between an ENM and a plant population's presence is propagule dispersal and establishment, not the distribution of the fundamental niche alone [@guisan2005predicting; @franklin2010moving; @engler2009migclim].
Previous innovations have largely focused on connecting currently suitable habitat to habitat predicted to be suitable under climate change [@engler2009migclim; @bocedi2014range], rather than on detecting populations that already exist but remain undiscovered.
To assist practitioners in detecting new populations or extending the range of known populations, we propose modelling plant occupancy as a random variable that depends on distance from known occupied sites.
'Distance' may be defined as Euclidean (or Haversine for large distances), or as a least-cost distance reflecting a generalised surface which conveys the difficulty for seeds to travel between the nearest occupied sites and the site of interest [@etherington2016least].
By combining habitat suitability with both Euclidean and least-cost distance, we distinguish candidate sites as extensions of known populations, new sub-populations, or new populations, rather than relying on suitability alone to prioritise field verification.

---

## P5 — Estimating population census size

- Obtaining reliable estimates of plant population census sizes requires two pieces of information: 1) measurements of plant density and 2) population extent.
- Often, instead of directly determining either the population's boundaries or density, observers estimate population sizes; however, these estimates can be unreliable [@ReischChristoph2018Acom].
- Several field methods for acquiring density measurements have been developed [@rominger2019application; @krening2021sampling; @ermakova2021densities; @alfaro2019optimal; @schorr2013using] to estimate population census size (*N*), as well as genomic methods to estimate effective population sizes (*N~e~*) (e.g. linkage disequilibrium) [@waples2024practical].
- However, these methods are focused on descriptive rather than predictive processes, estimating uncertainty in measurements *within* a population rather than across the species' range.
- We propose generating density estimates with spatial statistical models that predict both estimates and uncertainty across gridded surfaces [@oliver2012population; @doser2022spabundance; @a2022species].
- The second requirement, mapping a population's boundaries, is similarly time-consuming but no less essential for a realistic census size estimate.
- Population boundaries are currently delineated implicitly: practitioners walk outward until gene flow between individuals should be minimal (e.g., 1km), searching for additional individuals, or rely on aerial imagery [@alma9955306914202441].
- Many rare plant species are small and inconspicuous, and many taxa grow in small, clustered groups in specific habitats, making detection via either method difficult [@chen2013imperfect; @condit2000spatial].

### Draft prose

Obtaining reliable estimates of plant population census sizes requires two pieces of information: 1) measurements of plant density and 2) population extent.
Observers often estimate population sizes rather than directly measuring boundaries or density, and these estimates can be unreliable [@ReischChristoph2018Acom].
Several field methods for acquiring density measurements have been developed [@rominger2019application; @krening2021sampling; @ermakova2021densities; @alfaro2019optimal; @schorr2013using] to estimate population census size (*N*), as well as genomic methods to estimate effective population sizes (*N~e~*) (e.g. linkage disequilibrium) [@waples2024practical].
However, these methods are focused on descriptive rather than predictive processes, estimating uncertainty in measurements *within* a population rather than across the species' range.
We propose generating density estimates with spatial statistical models that predict both estimates and uncertainty across gridded surfaces [@oliver2012population; @doser2022spabundance; @a2022species].
The second requirement, mapping a population's boundaries, is similarly time-consuming but no less essential for a realistic census size estimate.
Population boundaries are currently delineated implicitly: practitioners walk outward until gene flow between individuals should be minimal (e.g., 1km), searching for additional individuals, or rely on aerial imagery [@alma9955306914202441].
Many rare plant species are small and inconspicuous, and many taxa grow in small, clustered groups in specific habitats, making detection via either method difficult [@chen2013imperfect; @condit2000spatial].

*Wordiness pass: "instead of directly determining either X or Y, observers estimate..." → "observers often estimate... rather than directly measuring X or Y".*

---

## P6 — All hypotheses together

- Here we: 1) evaluate the results of the best-performing model with groundtruthing (within 450m of existing occurrences).
- 2) Compare the results of models from different resolutions to the best model, as based on hold-out test data.
- 3) Determine whether model predictions are useful to find new sub-populations (450m - 1km from existing occurrences), or populations (> 1km from existing occurrences).
- 4) Estimate population census size using 4a) estimates of plant density, and 4b) modelling of population boundaries.

---

## Notes / gaps to fill during rewrite
- P2 (field validation) is currently thin on "across a species' range" — may need a new sentence about spatial coverage/representativeness of validation, not just ground-truthing volume.
- Several `[CITE]` / `(cite)` / `(taddeo2026diversitydistributions)` placeholders carried over — still need real citation keys before this becomes final text.
- The sentence "An ENM predicts a single variable..." and "Further, this progress... offers promise for modelling additional population parameters" both act as bridges — decide whether they open P1 (as a closer) or open P4/P5 as a topic sentence.


## orphaned


- An historic complication with the implementation and interpretation of ENMs is a mismatch between the spatial resolution of the independent variables available to model the species' fundamental niche and the factors governing the true distribution of populations - the realised niche [@carscadden2020niche; @chauvier2022resolution; @lembrechts2019incorporating].
 << orphan >>
- High-resolution spatial data eventually enable the determination of the extent of individual populations by detecting the edges of suitable habitat. *[also fits: P5 — population boundary mapping]* <<orphan>>

- A further mismatch of resolution is the year in which data on geographic localities were obtained and current conditions, which allow for positive population growth [@bracken2022maximizing].<<orphan>>
- Historic occurrence data may now reflect conditions inhospitable to population maintenance, or even represent populations that were simply sinks from more robust populations [@bracken2022maximizing]; using these records may decrease model performance.<<orphan>>
