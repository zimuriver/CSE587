(* Error messages *)
diceGibbs::oneFace = "The input to diceGibbs must contain at least 2 different faces.";
On[diceGibbs::oneFace];
 
 diceGibbs[input_, maxIterations_, outputFrequency_, outputFile_:"diceGibbsOutput.txt"]:=
	Module[{numFaces, binCounts, combinedBinCounts, thetaHyperParams, labelsStack={}, 
		    dirichletSample1, dirichletSample2, thetasStack={}},
		numFaces = Max[Flatten[input]];
		If[numFaces < 2,
		   Message[diceGibbs::oneFace];
		   Abort[]];
		(* Summarize each draw by the number of times each face occurred in it. *)
		binCounts = Map[BinCounts[#, {Range[numFaces+1]}] &, input];
		(* The total number of times each face occurred, summed over all draws.*)
		combinedBinCounts = Total[binCounts];
		(* The Dirichlet hyper parameters are the same for both die types. They are proportional
		   to the overall frequencies of the faces. Here we make them sum to 10, a modest prior
		   strength. *)
		thetaHyperParams = 10 * (combinedBinCounts / Total[combinedBinCounts]);
		(* labelStack is a list of all sampled label sets, from most recent to least recent. 
		   The initial labeling is random, with probability 0.5 for each label. *)
		labelsStack = {1 + RandomVariate[BernoulliDistribution[0.5], Length[input]]};
		dirichletSample1 = RandomVariate[DirichletDistribution[thetaHyperParams]];
		dirichletSample2 = RandomVariate[DirichletDistribution[thetaHyperParams]];
		(* thetaStack is a list of all the theta pairs sampled throughout the run, with the most recently 
		   sampled pair at the front. A pair is a list {theta1, theta2} where theta1 is a list of face 
		   probabilities for die type 1 and theta2 is a list of face probabilities for die type 2. 
		   The initial theta set is randomly sampled from the Dirichlet distribution based on the hyper
		   parameters. *)
		thetasStack = {{Append[dirichletSample1, 1-Total[dirichletSample1]],
					    Append[dirichletSample2, 1-Total[dirichletSample2]]}};
		(* Set up so that the output file appears in the same directory as the notebook from which this
		   function is called. *)
		SetDirectory[NotebookDirectory[]];
		Do[(
			PrependTo[labelsStack, 
			          sampleLabels[binCounts, labelsStack[[1]], thetasStack[[1]], {2, 2}]];
			PrependTo[thetasStack, 
				      sampleThetas[binCounts, labelsStack[[1]], thetaHyperParams]];
		    If[Mod[iteration, outputFrequency] == 0,
		       stackStats[labelsStack, thetasStack];
		       Put[labelsStack, thetasStack, outputFile]]),
		   {iteration, maxIterations}];
		{Mean[thetasStack], Mean[labelsStack]-1}
	]