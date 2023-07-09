# machine-learning
The machine learning framework consisted of ten separate algorithms and their combinations: random survival forest (RSF), elastic network (Enet), Lasso, Ridge, stepwise Cox, CoxBoost, partial least squares regression for Cox (plsRcox), supervised principal components (SuperPC), generalized boosted regression modeling (GBM), and survival support vector machine (survival-SVM).

Out of these, six algorithms (RSF, Enet, Lasso, Ridge, Stepwise Cox, and CoxBoost) included feature selection capabilities.

The RSF model was implemented using the randomForestSRC package, with parameters ntree (number of trees in the forest) and mtry (number of randomly selected variables for splitting at each node). The optimal values for mtry and nodesize parameters were obtained using the tune.rfsrc function and minimizing the out-of-sample error under 1000 ntrees.

The Enet, Lasso, and Ridge algorithms were implemented using the glmnet package. The regularization parameter λ was determined through 10-fold cross-validation, while the L1-L2 trade-off parameter α was set in the range of 0 to 1 (interval = 0.1).

The stepwise Cox model was implemented using the survival package, applying a stepwise algorithm based on the Akaike information criterion (AIC).

The CoxBoost model was implemented using the CoxBoost package, which fits a Cox proportional hazards model through componentwise likelihood-based boosting. The optimal penalty (amount of shrinkage) was determined using the 10-fold cross-validation routine optimCoxBoostPenalty. The number of boosting steps was selected using the cv.CoxBoost function, and the dimension of the selected multivariate Cox model was set using the principal routine in CoxBoost.

The plsRcox model was implemented using the plsRcox package. The cv.plsRcox function was used to determine the appropriate number of components, and the plsRcox function was applied to fit a partial least squares regression generalized linear model.

The SuperPC model was implemented using the superpc package, which is an extension of principal component analysis. It identifies linear combinations of the features that capture the largest variations in the dataset. The superpc.cv function employed a form of 10-fold cross-validation to estimate the optimal feature threshold in supervised principal components. To handle small validation datasets when fitting Cox models, it utilized the “pre-validation” approach.

The GBM model was implemented using the superpc package. The cv.gbm function employed 10-fold cross-validation to select the optimal number of trees with the minimum cross-validation error. The gbm function was then used to fit the generalized boosted regression model.

The survival-SVM model was implemented using the survivalsvm package. This regression approach accounts for censoring when formulating the inequality constraints of the support vector problem.

Each of these models used specific packages and methods to implement and optimize their respective algorithms.

# performance
Cox proportional hazards model and Kaplan-Meier analysis were performed with the ‘survival’ package. The receiver operating characteristic (ROC) curve was used to assess the prognosis classification performance of MALMS. The area under the curve (AUC) was calculated via ‘timeROC’ package. The comparisons between clinical and molecular traits, and risk score were implemented by ‘compareC’ package.

# mutation
Mutational landscape depiction and signatures extraction were both applied in the ‘maftools’ package. ExtractSignatures function based on Bayesian variant nonnegative matrix factorization factorized the mutation portrait matrix into two nonnegative matrices ‘signatures’ and‘contributions’, where ‘signatures’ represent mutational processes and ‘contributions’ represent the corresponding mutational activities[17]. The SignatureEnrichment function can automatically determine the optimal number of extracted mutational signatures and assign them to each sample based on the mutational activities. 

# drug sensitivity
The ‘oncoPredict’ package was used to build the drug sensitivity prediction procedure. The imputations were performed based on the expression matrix of a training set with known drug treatment information against the Genomics of Drug Sensitivity in Cancer (GDSC) database. The drug sensitivity scores of the samples were calculated using Ridge regression.

# MNF
The expression (Matrix A) was factorized into 2 nonnegative matrices W and H (i.e., A≈WH). Repeated factorization of matrix A was performed and its outputs were aggregated to obtain consensus clustering of samples. The optimal number of clusters was selected according to cophenetic, dispersion, and silhouette coefficients. The ‘NMF’ package with the brunet algorithm and 10 nruns algorithm was used to perform the consensus clustering.
