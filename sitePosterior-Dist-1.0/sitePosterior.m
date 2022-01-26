
(*A draw is the observed bases of a drawn sequence, which is subsequently replaced.
  sitePosterior calculates the posterior probability of a bound site versus a non-bound site, based on the
  observed sequence in the draw, the proportion of bound vs. non-bound sequences in the bag, and the 
  probabilities of observing each base in a sequence for bound and non-bound sites.
  The single number returned is the posterior probability of a bound site.*)
sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
Module[{trials,probOfSites,probOfBackground,sequenceCounts},
 
 	
 	probOfSites = Map[siteProbs[[#]][[sequence[[#]]]]&, Range[Length[sequence]]];
 	probOfBackground = Map[backgroundProbs[[sequence[[#]]]]&,Range[Length[sequence]]];
 	
 	
 	(Apply[Times,probOfSites]*sitePrior)/(Apply[Times,probOfSites]*sitePrior + Apply[Times,probOfBackground]*backgroundPrior)
    
    
		  
]
