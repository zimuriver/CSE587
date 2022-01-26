(* Wolfram Language package *)
(* Use this file to write code for supervised and unsupervised HMM parameter estimation.*)

(* supervisedHMMParameterEstimation should carry out maximum likelihood parameter estimation
   based on a set of observations and state labels, such as those provided in Test/mixed2.fa
   and Test/mixed2key.fa. It should output an HMM object that can pass checkHMMValidity. All
   the information needed to create an HMM object, including the number and names of the states
   and observation letters, should be taken from the two provided arguments. 
   
supervisedHMMParameterEstimation[observations_, stateSequence_]:=
*)

(* unsupervisedHMMParameterEstimation should carry out EM/ForwardBackward parameter estimation
   based on a set of observations. No state labels are provided. It should output an HMM object 
   that can pass checkHMMValidity. The number and names of the and observation letters, should 
   be taken from the first argument. The second argument is a list of strings that will serve
   as state names. It's length determines the number of states. 

unsupervisedHMMParameterEstimation[observations_, stateNames_]:=

*)