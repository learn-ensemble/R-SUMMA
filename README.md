# SUMMA (Spectral Unsupervised Multiple Method Aggregator)
SUMMA is an R package for unsupervised ensemble learning which combines predictions of multiple predictors.
Ensemble learning is an elegant solution for the over-fitting problem in machine learning which is an important
challenge in computational biology. The biggest advantage of SUMMA is that it is unsupervised so it does not require gold standard
labels associated with each sample which is not available in most biomedical problems.
Therefore, SUMMA can be applied to a wide range of biomedical problems including Network Inference,
Differential Expression, Cancer Diagnosis or any other two class classification problem for which there is
more than a single algorithm is available.

SUMMA needs as an input a prediction matrix, denoted as  *P*, which represents the predictions
of individual methods.
The matrix *P* is of size *m* by *n*, where *m* is the number of samples (observations) and *n* is the number of methods. SUMMA supports two modes for the prediction matrix it can either have binary numbers *{-1,1}* in which case SUMMA will run the SML algorithm of [Parisi et.al.](http://www.pnas.org/content/111/4/1253.abstract),
or *P* may consists of real numbers that allows the algorithm to rank samples, in which case the SUMMA algorithm will be run on the data. In the rank prediction case by convention we assume the positive samples are more likely to have higher scores than samples belonging to the negative class. With the prediction matrix as an input, SUMMA will do the following:

1. SUMMA will rank methods based on their performance:
    - If the matrix *P* has binary values, the ranking will be based on SUMMA estimated balanced accuracies of the individual methods
    - If the matrix *P* has continuous values, the ranking of the methods will be based on SUMMA estimated AUCROC of individual methods

2. SUMMA will calculate weights for each methods proportional to its performance which also corresponds to the maximum likelihood estimator. Using the weights SUMMA will calculate an unsupervised ensemble classifier.

For more details about the theory of the SUMMA algorithm please check our paper [Ahsen et. al.] (https://arxiv.org/abs/1802.04684).
In the SUMMA package, we also include the so-called WOC (Wisdom of Crowds) ensemble which assigns equal weight to each method.
The WOC ensemble has shown robust performance throughout various crowdsourcing challenges such as the DREAM5 Network Inference
Challenge [Marbach et. al.] (https://www.nature.com/nmeth/journal/v9/n8/abs/nmeth.2016.html).
Next two sections contains some example codes for the binary and rank prediction cases.

## Binary predictions
For the convenience of the user, we have added `create_predictions(m,n,p,type)`, where *m* is the number samples, *n* is the number of methods, *p* is the prevalence and *type* is type of output which either
the ``"binary"`` or `"rank"`. We can create a binary classification problem with `m=3000,n=30` and
`p=0.3`.
```r
data_binary=create_predictions(3000,30,0.3,"binary")
```
For consistency issues, we will load the
data included with the package
```r
load(system.file("extdata", "data_binary.Rdata", package="summa"))
prediction_matrix=data_binary$predictions
actual_labels=data_binary$actual_labels
```
The function `create_predictions` outputs a list where the first element of the list is the prediction matrix and second element is the gold standard labels associated with each sample. Note that the actual labels is not necessary to run SUMMA, however if you have it SUMMA has a convenience function to evaluate the results for the user. Next let us run summa:
```r
summa=summa(prediction_matrix,"binary")
```
which gives an S4 class summa object.
The most relevant parts of the summa object is `summa@estimated_performance` which is the estimated balanced accuracy, and `summa@estimated_label` which is the estimated labels for each sample. If we also have the gold standard labels we can use the `calculate_performance` function to calculate the actual performances as well as the performance of the unsupervised ensemble classifier. To do this we run
```r
summa=calculate_performance(summa,actual_labels)
```
which returns again a summa object in which `summa@actual_performance` performance is the actual performance of each method and `summa@summa_performance` is the performance of the unsupervised ensemble. Moreover, `summa@woc_performance` corresponds to the performance of the WOC (wisdom of the crowds) where one assign equal weights to each classifier. Finally, we have also added a convenience plot function which plots the actual v.s. estimated performances of each method as well as the performance of unsupervised ensemble and WOC
```r
summa_plot(summa)
```

![binary_erf](/figures/performance_binary.png)
## Ranked Predictions
When the prediction matrix have continuous values that allows the samples to be ranked, the package runs the summa algorithm. The intuition and outputs in this case
are exactly the same as in the binary case but in this case the performance is the AUCROC instead of balanced accuracy. Again we create predictions using `create_predictions` functions but this time we set
`type=''rank''` and also we can create a prediction
matrix with `m=2000` samples and `n=40` methods with
prevalence  `p=0.7`
```r
data_rank=create_predictions(2000,40,0.7,"rank")
```
We will load the
data included with the package
```r
load(system.file("extdata", "data_rank.Rdata", package="summa"))
prediction_matrix=data_rank$predictions
actual_labels=data_rank$actual_labels
```
Next we run SUMMA using the simulated data
```r
summa=summa(prediction_matrix,"rank")
```
where the only change from the binary case is that we set the second argument as `"rank"`.
If we have the actual labels, we can again use the
`calculate_performance` function to calculate the performances
```r
summa=calculate_performance(summa,actual_labels)
```
where again `summa@estimated_performance`  is the estimated AUCROC, and `summa@estimated_label` is the estimated labels for each sample. Different from the binary case, `summa@estimated_rank` gives a ranking of the samples predicted by the summa algorithm. This might be useful for comparing different samples. Finally, we can use the `summa_plot` function to plot the performances
```r
summa_plot(summa)
```
![binary_erf](/figures/performance_rank.png)
