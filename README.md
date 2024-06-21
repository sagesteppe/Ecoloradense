## Iterative Adaptive Modelling - a test with Eriogonum coloradense

Eriogonum coloradense has a bimodal distribution, with most populations occuring in alpine habitats on a dozen mountains in the central Southern Rocky mountains of Colorado, with a handful of small populations occuring in sagebrush meadows to the south and East. 

Species Distribution Models were made for Eriogonum colardense using existing occurrence records obtained via herbaria consortia (Southern Rockies Herbarium Consortium), and iNaturalist. 
While iNaturalist records are ostensibly associated with relatively accurate (+/- 10m) GPS coordinates, herbarium records are generally less accurately mapped. 
The first iteration SDM's captured the coarse ecological niche (Figure X) of the species, but lacked resolution distiguishing the species more generally occupied habitat on rocky, exposed, ridgelines. 
These models were used to guide field sampling efforts that had two efforts, verify the accuracy of all existing occurrence data as encountered, and two generate a much larger data set of occurrences and localized absences for the second modelling iteration. 

Spatial data products of all trails known to the analysts in the alpine population cluster, and found via review of topographic maps and satellite imagery, were generated to define areas which were feasible for field sampling. 
500 random points, all separated by more than 90meters, within 100 meters of these locations were then generated for ground truthing. 
Given the noted relative sparsity of the taxon any occurrences encountered had data collected at their sites. 

## Methods

All herbarium records used as occurrences for E. coloradense were manually cleaned by viewing the geolocation along with label data, records with discrepancies were discarded. 
Records were then thinned to remove duplicates, and then thinned again using spThin with a mininum distance of 90m, and one final subset containing the most records, 81, was randomly selected. 
