
(*A draw is the series of rolls between the time a die is drawn from the bag and the time it is returned to the bag.
  dicePosterior calculates the posterior probability of a type1 versus type2 die, based the number of times each face
  appears in the draw and the relative numbers of Type 1 and Type 2 dice in the bag as well as the face probabilities
  for Type 1 and Type 2 dice. The single number returned is the posterior probability of Type 1.
  
  You will need to catch the case where the probability of a face showing is zero and it shows zero times. 0^0 is undefined,
  but there is a simple intuitive number that you should use in this case.
  *)

dicePosterior[binCounts_, type1Prior_, type2Prior_, faceProbs1_, faceProbs2_] := Module[{trials,probabilityResultsGiven1,probabilityResultsGiven2},
 trials = Total[binCounts];
 
 If[Length[binCounts]==Length[faceProbs1]==Length[faceProbs2],
 	probabilityResultsGiven1 = PDF[MultinomialDistribution[trials, faceProbs1], binCounts];
 	probabilityResultsGiven2 = PDF[MultinomialDistribution[trials, faceProbs2], binCounts];

 	(probabilityResultsGiven1*type1Prior)/(probabilityResultsGiven1*type1Prior + probabilityResultsGiven2*type2Prior),
    
     Return["Bin counts and face probabilities must match"]] (*If improper input is given*)
 ]
