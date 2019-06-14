# Simulation pipeline instructions	

Here we provide a customizable tool for simulation of evolve and resequence experiments on a quantitative trait to evaluate their power for QTL detection. The simulation is carried out in the program SLiM and the data analysis is conducted in R.

Below is our suggested workflow. Each step in the workflow has to be completed before the next can start.

For the simulation of E&R experiment:

1.	Create neutral populations representing the population prior to the selection experiment. To do this, read and edit `ShellScripts/Burnin.sh`. Then, run `Burnin.sh`.

2.	Establish the trait architecture and simulate the selection experiment. For this, read and edit `ShellScripts/Selection.sh`, then run the script. This step can be repeated to simulate different combinations of trait architecture and experimental design. These can be created using the same burn-in.

For the power analysis in QTL detection:

3.	Compile the SLiM outputs from each simulation in step 2, and (optionally) create input files for WFABC and ApproxWF. To do this, read, edit, and run `RScripts/Compile.R`

4.	Calculate the power and false positive rate for each simulation. For this, read, edit, and run `RScripts/Analysis.R`

5.	Plot ROC curve comparisons by combining the ROC tables from simulations of different trait architecture or experimental design given by step 4. You can use `RScripts/PlotROC.R` to do this. Alternatively, you can also create your own R scripts, since this step highly depends on what your specific comparison is.

Note: 

To adjust most of the variables, you will not need to edit the SLiM script and can directly input them through `Burnin.sh` and `Selection.sh`. Below is a list of these variables.

|Trait architecture	|Population parameters	|Experimental design setting|
|---|---|---|
|Number of QTLs	|Population size*	|Sample size|
|Position of QTLs	|Number of chromosomes	|Length of experiment*|
|Effect sizes of QTLs	|Length of chromosomes	|Mode of selection|
|Starting frequency of QTLs	|Recombination rate	|Strength and direction of selection|
|Dominance	|Nucleotide diversity	| |
|Pairwise epistasis	|	| |

\* To modify these two variables, one single value will need to be changed in the `.slim` files. See the scripts for details. 

There are some variables and scenarios that are not explicitly incorporated in the SLiM script but that can be conveniently implemented by editing `SlimScripts/Burnin.slim` and `SlimScripts /Selection.slim`. These variables and scenarios include:

* heritability values other than 1 (see the SLiM manual 13.4)
* selection modes other than truncating selection (see the SliM manual Chapter 13.1)
* pleiotropy (see the SLiM manual Chapter 13.5)
* population structure (see the SLiM manual Chapter 5.2)
* demographic history prior to the selection experiment (see the SLiM manual Chapter 5)
* known genotypic data of the experimental population (see the SLiM manual Chapter 18.12)

For data analyses and visualizations (step 3-5), it might be more efficient if you write your own scripts from scratch. 

If WFABC and/or ApproxWF are used, they should to be run between step 3 and 4. 

Please post an issue on GitHub or contact us at rl683@cornell.edu with any problems or questions.
