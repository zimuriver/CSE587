(* makeCVFolds calls makeCVTestFolds which divides the list up into numFolds disjount subsets whose lengths differ
   by at most one. Each subset, or "fold" is used as a held-out test data for a different run. For each test subset, 
   makeCVFolds makes a pair consisting of the training fold (everythign that's not in the test fold) and the test fold. 
   A list of such pairs, one per test fold, is returned.
   
   The following example is the result of running makeCVFolds on {a, b, c, d} and 4.
   
   {
    {{b, c, d}, a},
    {{c, d, a}, b}
    {{d, a, b}, c}
    {{a, b, c}, d}
    }
   *)
   
makeCVFolds[list_, numFolds_]:= (
   If[Length[list] < numFolds,
   	  Print["makeCVFolds: The number of folds, ", numFolds, " cannot exceed the length of the list, ", Length[list]];
  	  Print["List: ", list];
   	  Abort[]];
   Map[{Complement[list, #], #} &,
   	   makeCVTestFolds[list, numFolds]]
   	   )

(* makeCVTestFolds which divides the list up into numFolds disjount subsets whose lengths differ
   by at most one.*)

makeCVTestFolds[list_, numFolds_]:=
 	With[{length = Length[list]},
 	  TakeList[list, Catenate[{ConstantArray[Floor[length / numFolds], 
 		                                     numFolds - Mod[length, numFolds]],
 		                       ConstantArray[Floor[length / numFolds] + 1, 
 		                                     Mod[length, numFolds]]		                     
 		                     }]]]

(* constructEdgeList takes an edgeMatrix with one row (inner list) for each target gene and one 
   column for each predictor variable. The entries are the scores from the final regression. 
   ConstructEdgeList outputs a list of triplets, {predictorIndex, targetIndex, scores}, one 
   triplet for each non-zero score in the input matrix. *)

constructEdgeList[edgeMatrix_] :=
   Module[{predictorIndex, targetIndex,
   	       targetCount=Dimensions[edgeMatrix][[1]],
   	       predictorCount=Dimensions[edgeMatrix][[2]],
   	       edgeList={}},
   	  For[targetIndex = 1, targetIndex <= targetCount, targetIndex++,
         For[predictorIndex = 1, predictorIndex <= predictorCount, predictorIndex++,
         	If[Abs[edgeMatrix[[targetIndex, predictorIndex]]] > 0,
         	  AppendTo[edgeList, 
         	  	       {predictorIndex, targetIndex, edgeMatrix[[targetIndex, predictorIndex]]}]
         	  ] 
         	]
         ];
      edgeList]
 
 (* sortCSVRowColumn makes sure all the matrices have their genes in the same order. 
    Note that strings including letters and numbers will be sorted in dictionary
    order, going one character at a time through the string. Thus, if strings
    contain numbers they will not necessarily be sorted in numeric order.*)

sortCSVRowColumn[labeledMatrix_] := 
   Module[{sortedRows, sortedTransposeRows},
   	 (* Sort the rows by row labels, except for the first row which consists
   	    of column labels. *)
   	 sortedRows=SortBy[labeledMatrix[[2 ;;]], First];
   	 (* Add the column labels back*)
   	 PrependTo[sortedRows, labeledMatrix[[1, All]]];
   	 sortedTransposeRows = SortBy[Transpose[sortedRows][[2 ;;]], First];
   	 Transpose[Prepend[sortedTransposeRows, sortedRows[[All, 1]]]]
   ]
  
(*myRound is a hack to get around a problem with rounding numbers in Mathematica.Please don't bother to try to understand it.*)
myRound[x_, n_] :=
       N[IntegerPart[Round[x, 10^-n]*10^n]
                            /10^n];	                     

(* There is a separate regression to predict the expression of each gene, but the predictors for each gene are the same.
   The only difference is that, when predicting the expression of a TF gene, you don't include that gene as a predictor.
   However, that difference is handled elsewhere. constructPredictorsMatrix constructs a matrix whose columns are the 
   expression levels of all the TFs and whose rows are the samples. This is the first step to constructing the regression
   problem.
   
   In the predictors matrix (design matrix) each row is a sample and each column is a TF. The entry is the expression 
   level of the TF in the sample. The expression Matrix must be transposed and the entries corresponding to TF genes
   selected out. *)
