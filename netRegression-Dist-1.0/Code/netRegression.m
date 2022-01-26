(* Implements TF network inference using the LASSO, a.k.a. L1-penalized regression approach. *)

(*
Basic outline of tasks:
-Read in the target expression matrix (genes x samples)
-Read in the TF expression matrix (TFs x samples)
-Construct folds for cross-validation
-Map over all the target genes
---map over lambdas
-----map over folds Fit[{predictors, responses}, FitRegularization -> {"L1", lambda}]
---calculate mean error across folds for the current lamda
---Choose lambda
---Run Fit again on the full data set
---Extract edge scores 
*)

(* inferTFNetLasso does all the stuff that only needs to be done once for all target genes:
-Import the predictor (TF) expression matrix from a .csv into a matrix.
-Import the target gene expression matrix from a .csv into a matrix.
-Construct predictor folds for cross-validation (these are the same for all genes)
-Map over all the target genes, passing in the information just computed plus the
 CV folds of one target gene.
*) 

inferTFNetLasso[predictorsMatrixPath_, geneExpressionMatrixPath_, cvFolds_:10, lambdas_:{0}] :=
	Module[{netObject, predictorFolds,
		    predictorsInput=sortCSVRowColumn[Import[predictorsMatrixPath]],
		    geneExpressionInput=sortCSVRowColumn[Import[geneExpressionMatrixPath]]},
	(* check that the sample labels in the two input files are the same.*)
	If[predictorsInput[[1, 2;;]]!=geneExpressionInput[[1, 2 ;;]],
	   Print["The sample labels in the gene expression matrix do not match those in the predictors matrix."];
	   Abort[],
	   netObject["samplesNames"]=predictorsInput[[1, 2;;]]];
	(* Transpose so that samples (training examples) are rows and features (TFs) are columns. *)  
	netObject["predictorsExpressionMatrix"]=Transpose[predictorsInput[[2 ;;, 2 ;;]]];
	netObject["predictorNames"]=predictorsInput[[2 ;;, 1]];
	(*Do not transpose -- leave genes as rows, samples as columns. Each row will be used in a separate regression. *)
	netObject["geneExpressionMatrix"]=geneExpressionInput[[2 ;;, 2 ;;]];
	netObject["geneNames"]=geneExpressionInput[[2 ;;, 1]];
	(* Results of LASSO regression for each target gene will be stored in netObject["edgeScores"]. *)
	netObject["edgeScores"]={};
	netObject["lambdasChosen"]={};
	(* makeCVFolds is defined in tools.m Read the description to understand what it returns. *) 
	predictorFolds = makeCVFolds[netObject["predictorsExpressionMatrix"], cvFolds];
	(* netLasso actually does the regression for one target gene and returns the vector of regression coefficients.*)
	Map[netLasso[netObject, predictorFolds, makeCVFolds[#, cvFolds], lambdas] &,
	    netObject["geneExpressionMatrix"]];
	(* The matrix of coefficients is now converted to an adjacency list containing (TF, target) pairs that have
	   non-zero coefficients, sorted by the absolute value of the coefficient from largest to smallest.*)
	netObject["rankedEdges"] = ReverseSortBy[constructEdgeList[netObject["edgeScores"]], Abs[Part[#, 3]] &];
	netObject	
	]
   
 (* netLasso takes in a netObject, where inputs and results are stored, predictor values (TF expression 
    levels in each sample) and response values (target gene expression levels of one gene in each sample). 
    The predictor and response values are divided into folds -- see the documentation on makeCVFolds in tools.m. 
    The last input is a list of lambda values to try.
 
    netLasso calls constructErrorsMatrix and passes the result to chooseLambda, which returns the lambda value
    that minimizes cross-validation error.
    
    Once the best lambda value has been chosen, netLasso calls the built in function Fit on the complete predictors
    matrix, including all samples, the complete response vector, also including all samples, and the best lambda.
    "Fit" does the regression with the chosen lambda value and returns the coefficients. Please read the docs on
    Fit and note that the term 'design matrix' means the same thing as 'predictors matrix', in our case a matrix
    whose rows are samples and columns are TFs. The resulting coefficients are the edge scores for all the TFs on 
    this gene and are added to the end of the matrix of edge scores as a final row.
    *)

netLasso[netObject_, predictorFolds_, geneFolds_, lambdas_]:= (

       AppendTo[netObject["lambdasChosen"],
       chooseLambda[constructErrorsMatrix[predictorFolds,geneFolds,lambdas],lambdas]
       ];

       (* Fill in here code that returns the best lambda. To do that, you'll need to call constructErrorsMatrix
          and call chooseLambda on the result. *) 

       AppendTo[netObject["edgeScores"],
       Normal[Fit[{Join[predictorFolds[[1,1]],predictorFolds[[1,2]]], Join[geneFolds[[1,1]],geneFolds[[1,2]]]}, "BestFitParameters", FitRegularization->{"L1", netObject["lambdasChosen"][[-1]]}]]
       ];
      (* Fill in here code that returns the coefficients of the best fit, given the best lambda. To do that, you'll 
         need to call Fit with the complete predictor set (no folds), the complete response set (no folds), and the 
         best lambda. You can pull the best lambda back from the netObject with netObject["lambdasChosen"][[-1]]. *) 

      (* Return the modified netObject. *)
     
      netObject
      )

(* The errors matrix has one row (inner list) for each possible value of lamda. Each row has one entry for each
   fold used in cross-validation. The entry is the sum of squared residuals when using that lambda value on that fold. *)

constructErrorsMatrix[predictorFolds_, geneFolds_, lambdas_]:= 
    (* You must map or iterate over lambdas and then folds. For each (lambda, fold) pair you will call sumSquaredTestError,
       which calls Fit to do the actual penalized regression on a given training fold and then calculates the prediction
       error of the model on the test fold. *)
	Module[{errorsMatrix,i,j},
	

	errorsMatrix = ConstantArray[0,{Length[lambdas],Length[geneFolds]}];
	Do[
		Do[
			errorsMatrix[[i,j]] = sumSquaredTestError[predictorFolds[[j]],geneFolds[[j]],lambdas[[i]]],
				{i,Length[lambdas]}],
					{j,Length[geneFolds]}];
	
	errorsMatrix
	
	]

(* sumSquaredTestError calls Fit to do the actual penalized regression on a given training fold and then calculates the 
       prediction error of the model on the test fold. *)
sumSquaredTestError[predictorFold_, responseFold_, lambda_] :=
   Module[{predictorTrain=predictorFold[[1]],
   	     predictorTest=predictorFold[[2]],
         responseTrain=responseFold[[1]],
   	     responseTest=responseFold[[2]],
   	     bestFitParameters 
   	    },
   	  (* THIS IS WHERE THE REGRESSION FOR OPTIMIZING LAMBDA ACTUALLY HAPPENS. *)
   	  bestFitParameters= Fit[{predictorTrain, responseTrain}, "BestFitParameters", FitRegularization->{"L1", lambda}];
   	
      Total[(Map[Dot[bestFitParameters, #] &, predictorTest] - responseTest)^2]
   ] 

 
(* chooseLambda takes an errors matrix and computes the total error across folds for each lambda and then choses the 
   lambda whose error is least.
*)
chooseLambda[errorsMatrixLambdasByFolds_, lambdas_]:=
    (* Put your code for chooseLambda here. *)
   Module[{meanErrors,chosenLambda,i},
   		meanErrors = ConstantArray[0,Length[errorsMatrixLambdasByFolds]];
  		Do[meanErrors[[i]] = Mean[errorsMatrixLambdasByFolds[[i]]],{i,Length[errorsMatrixLambdasByFolds]}];

   		chosenLambda = lambdas[[Ordering[meanErrors,1][[1]]]];
   		
   	chosenLambda
   ]
(* Creates a precision-vs-rank plot and a precision-recall plot, with a line indicating
   the expected precision for a random selection of possible edges.
   which is the number of binding edges divided by the number of possible edges.*)
presentNetAccuracy[netObject_, bindingMatrixPath_] :=
	Module[{rankPrecisionRecall = scoreNet[netObject, bindingMatrixPath],
		    bindingImport = sortCSVRowColumn[Import[bindingMatrixPath]],
		    bindingMatrix, bindingEdgeCount, possibleEdgeCount, predictedEdgeCount, 
		    expectedRandomPrecision,
		    precisionRankPlot, precisionRecallPlot},
      bindingMatrix = bindingImport[[2 ;;, 2 ;;]];
      bindingEdgeCount = Total[bindingMatrix, 2];
      possibleEdgeCount = Apply[Times, Dimensions[bindingMatrix]];
      predictedEdgeCount =  Length[rankPrecisionRecall];
      expectedRandomPrecision = N[bindingEdgeCount / possibleEdgeCount];
      precisionRankPlot = 
      	Labeled[Show[ListLinePlot[rankPrecisionRecall[[All, 1 ;; 2]], PlotRange -> All], 
                     Graphics[Line[{{0, expectedRandomPrecision}, {predictedEdgeCount, expectedRandomPrecision}}]],
                     Graphics[Style[Text["Random Expectation", {predictedEdgeCount * 0.67, expectedRandomPrecision + 0.035}], Medium]]], 
                {Text["Edge Rank Threshold"], Rotate[Text["Precision (Fraction of edges supported)"], Pi/2]}, 
                	{Bottom, Left}];
      Print[precisionRankPlot];
      precisionRecallPlot =
      	Labeled[Show[ListLinePlot[Reverse[rankPrecisionRecall[[All, 2 ;; 3]], {2}], PlotRange -> All], 
                     Graphics[Line[{{0, expectedRandomPrecision}, {predictedEdgeCount, expectedRandomPrecision}}]],
                     Graphics[Style[Text["Random Expectation", {rankPrecisionRecall[[-1,-1]] * 0.67, expectedRandomPrecision + 0.035}], Medium]]], 
                {Text["Recall"], Rotate[Text["Precision (Fraction of edges supported)"], Pi/2]}, 
                	{Bottom, Left}];      
      Print[precisionRecallPlot];
	] 

scoreNet[netObject_, bindingMatrixPath_] :=
   Module[{bindingImport = sortCSVRowColumn[Import[bindingMatrixPath]], 
   	       bindingMatrix, targetGeneNames, TFNames},
      bindingMatrix = bindingImport[[2 ;;, 2 ;;]];
      targetGeneNames = bindingImport[[2 ;;, 1]];
      TFNames = bindingImport[[1, 2 ;;]];
      If[targetGeneNames != netObject["geneNames"],
      	 Print["The target gene names in the gene expression matrix ",  netObject["geneNames"], " do not match those in the binding matrix, " targetGeneNames];
      	 Abort[]];
      If[TFNames != netObject["predictorNames"],
      	 Print["The TF names in the TF expression matrix do not match those in the binding matrix."];
      	 Abort[]];
      evalEdges[netObject["rankedEdges"], bindingMatrix]
   ]
 
 (* evalEdges takes the list of ranked edges with non-zero scores that constructEdgeList outputs (see tools.m).
    Each edge has the form {TF#, target#, score} and they are sorted by the absolute value of the score, from highest to lowest. A
    positive score indicates that the TF activates transcription of the target gene, and a negative score indicates that the TF represses
    transcription of the target gene. The absolute value of the score indicates the "confidence level" or "effect size", which I 
    put in quotes because they are intuitive ideas of how we should use the score rather than mathematical definitions. (TF, target)
    pairs with score zero, which should be the vast majority, are not included in the ranked edges.
    
    evalEdges returns a different list of triplets, one for each of the input edges: {rank, precision, recall}. The rank is an integer 
    that indicates the position of the edge in the rankedEdges list. It is not strictly necessary, since the triplets are still in rank
    order, but it's convenient to have it in the triplet. The other two are precision and recall when the rank of the triplet is used
    as the threshold between the edges that are predicted to be real and those whose score is too small to be considered real. In other 
    words, the rank of the triplet is used to effectively binarize the absolute values of edge scores into one for "predicted" or 
    zero for "not predicted".
    
    Before you start, make sure you understand how precision and recall are calculated. To calculate them, you will need to know how
    many of the set of possible edges are supported by the binding data -- in  other words, how many ones are present in bindingMatrix.
    You should only calculate this once, for all ranks. You will also need to calculate the number of "true positives" above the rank
    threshold -- the number of higher ranked edges that are supported by ones in the binding matrix. I strongly suggest that you calculate
    the number of true positives by using a For loop to iterate through the ranks, starting at 1 and going down. Keep a running count
    of the number of true positives and increment it whenever you pass an edge that is supported by the bindingMatrix.
    
    The bindingMatrix input is a binary matrix whose rows are target genes and columns are TFs. Each entry is a one if the binding data
    support an edge between the corresponding TF and target; otherwise it is zero.
    
    Note that the order of TFs and target genes in bindingMatrix is guaranteed to be the same as in the TF expression matrix and the 
    target gene expression matrix. So if an edge in the rankedEges input is {7, 3, X} (TF number 7, target gene number 3, score X), 
    then the information about whether the edge is supported by the binding data is located in row 3, column 7, of bindingMatrix. 
    Note that the index order is switched from the edge to bindingMatrix.      
    *)     	 
evalEdges[rankedEdges_, bindingMatrix_]:=
   (* Put your code for evalEdges here.*)
 Module[{truePositives=0,tripletOutput,i,precision,recall,numberOfRankedEdges=Length[rankedEdges],numberOfRelevant},
tripletOutput = ConstantArray[{0,0,0},numberOfRankedEdges];
numberOfRelevant = Total[bindingMatrix,Infinity];
 	
	Do[
		If[bindingMatrix[[rankedEdges[[i,2]],rankedEdges[[i,1]]]]==1,truePositives++];
		precision = N[truePositives/i];
 		recall = N[truePositives/numberOfRelevant,MachinePrecision]; 	
		tripletOutput[[i]] = {N[i,MachinePrecision],precision,recall}
	,{i,numberOfRankedEdges}];
 
 tripletOutput
 ];