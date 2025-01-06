## Iterative Adaptive Modelling - a test with Eriogonum coloradense

Goals:

1) Evaluate the differences between species distribution models of 4 different resolutions (1m, 3m, 1/3, 1, and 3 arc-second) on both held-out and ground verified data [AUC-ROC, AUC-PR, Sens. Spec. TSS].  
2) Estimate spatial plant counts across species domain, using XGBoost, ADABoost and Poisson GLM. 
3) Estimate population census sizes using a combination of thresholding of SDM results, and count predictions. 
4) Model occupancy of suitable habitat patches as a function of distance from an occupied patch (both euclidean and least-cost), patch (e.g. inner area), and class (e.g. isolation) characteristics. [multiple logistic regression]. 
5) Compare the results of a Mature plant presence versus a Juvenile plant presence model, are the important  variables similar? Are patch sizes similar? 
6) Simulate the effects of sample size on model performance using the data set generated after both rounds of field sampling. 
7) Simulate the effects of occurrence record geolocation accuracy on model performance using the data set generated after both rounds of field sampling. 