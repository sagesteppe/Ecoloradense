## Iterative Adaptive Modelling - a test with Eriogonum coloradense

Goals:

1) Evaluate the differences between models of 5 different resolutions (1m, 3m, 1/3, 1, and 3 arc-second) on both held-out and ground verified data [AUC-ROC, AUC-PR, Sens. Spec. TSS].
2) Compare Results from 150 ground truth points along the axis of suitability (0-100%), and species occupancy (0-100% percentile of cells). Analyze as two independent logistic regressions, and together as two multiple logistic regressions.  
3) Model occupancy of suitable habitat patches as a function of distance from an occupied patch (both euclidean and least-cost), patch (e.g. inner area), and class (e.g. isolation) characteristics. [multiple logistic regression]. 
4) Compare the results of a Mature plant presence versus a Juvenile plant presence model, are the important  variables similar? Are patch sizes similar? 
5) Simulate the effects of sample size on model performance using the data set generated after both rounds of field sampling. 
6) Simulate the effects of occurrence record geolocation accuracy on model performance using the data set generated after both rounds of field sampling. 