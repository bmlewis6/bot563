2/7
Needleman-Wunsch
F(i,j)=min{F(i−1,j−1)+cost(ai,bj),F(i−1,j)+c,F(i,j−1)+c}
c= cost of gap

Align S1= ATCG, S2= TCA
  _ A T C G
_ 0 2 4 6 8
T 2 1 2 4 
C 4
A 6
(0->1), (2->2), (2->4)
arrow points at gap when vert/ho

Prgressive alignment relies heavily on guide tree, cannot fix alignment from leaves to roots (first alignment stays)
Profile- fraction of frequency of neucleotide at the position

T-Coffee Notes
Slower but less greedy bc it considers all alignments against each other.
Local similarity- when two proteins share only a domain or motif
Progressive alignment like ClustalW
The primary library can have multiple versions of aligning two sequences. All alignments are weighted
Library extensions, residues, triplets
position specific storing scheme

ClustalW downloading
download clustalw-2.1-win.msi and install
move to usalable folder (like bot563)
open the .exe program
put data in clustal folder
load data using option1
align with option2
then align with option2
name output file and guide tree file

2/21
character based- keep information on neulceotides
p-distance- fraction of nucleotides which are dissimilar between two sequences
possible to have multiple mutations hidden by last state (underestimating) -> apply evolutionary model
don't need to search space of trees- distance based returns the optimum tree
ultrametric- distance from root is the same to all tips
molecular clock- genetic distance is linear with time (rate of evolution is constant)
minimum evolution uses nj (same distance input as UPGMA)
neighbor joining corrects for unequal rates
N= number seq
can get neg rates in rate corrected distance matrix

2/23
Parsimony is character based, minimizes amount of evolutionary change, assumes rate is slow
(((A,B),(C,D)), E)
length of ((G,A),(C,C))- 2 by AAC configuration
dynamic programming- breaking a large task into smaller ones
parsimony score for (((C,C),T),G) {C} {C,T} {C,T,G} 2
if there is an imbalance in branch length, tree will come out wrong bc it will group shorter, known branches together

3/2
P(X=k)= (h^k * e^-h)/k!   h=lambda=avg k=number wanting prob of
avg= 8.4, P(X=5)=.07837 P(X=0)=.00022
assumptions and quality of date are most important
assumption 1- mutation process is the same at every branch-> look at one mutation to determine pattern
assumption 2 - sites evolve independently -> focus on mutation between two sites
assumption 3- all sites evolve the same-> choose any site to model the mutation
ut instead of lambda- rate of mutation events(mu)*expected number of mutation events in time t (mu t)
probability matrix is sum over all possible numbers of k = e^Qut Q being instantaneous rate matrix generator (R-I) I being identity matrix
homogeneous continuous-time markov chain
homogeneous- prob for next event doesn't depend on time point where chain is
continuous-time- event can happen at any time point (not just discrete steps)
Markov- prob of the next state only depends on current state
choosing substitution model means choosing Q
choose the simplest model possible bc fewer parameters is more accurate

3/7
substitution- a mutation to a nucleotide other than the orginal (A->A and A->G is mutation, A->G is substitution)
pi matrix- frequency of nucleotide occurance plus 0s
trace- sum of diagonal
Q and mu are restircted variables so the branches are scaled properly
among-site rate variation- sites have dif rates, assumes rate is constant
site-specific rate variation- site has a rate which can change later

3/9
for phylogenetic inference choose distance, parsimony, likelihood then search space for optimal tree
likelihood depends on model, branch length, and ancestral states
felsensteins pruning algorithm- have prob of observiting each nulceotide/aa at each node in tree, formula for left and formula for right, x interates on value of child, i iterates on value of parent, Lsub(l) is length of branch leading to left
RAxML starts with most parsimonious tree
factors affecting ML performance: starting tree, model chosen, data, convergence

3/28
Bayesian inference
MLE (max lik est) for lambda- need to maximize likelihood function by averaging sum of data over days
prior + likelihood = posterior
prior will pull posterior in direction of prior's mean
as you acquire more data, posterior migrates towards ML
likelihood x prior is proportional to posterior
conjugate prior- end up with same dist, dif parameters
mcmc- start at random lambda naught, for i propose new lambda (iterate), keep track of all values of lambda (chain), chain stays in lambda of high value posterior so histogram is approximation of posterior
mcmc considerations- choice of priors (use priors other have used), mixing (want chain to traverse whole space), convergence(how well do you reach regions of high posterior), burnin (how long does it take to reach these regions/region before high post)
smaller variance means takes longer to traverse space, but almost every step is accepted
white noise in trace plots is a good sign, asses whether chain has converged or not, long periods of time stuck  in one place indicates too large of steps
hot chain cold chain robots

RevBayes can be used to create model if existing method doesn;t work for your data

4/11 Coalescence model
at least 1 mil gen
internal branches of true tree are v short can lead to concatination converging wrong
cannot have recombination
gene tree --(coalescence)--> species tree
coalescence in one pop- assume no selection and one parent, uniform prob (1/N)
coalescent event- when two ind have common ancestor 
P(coalescence)= 1/N
P(no coalescence)=1- 1/N
P(no coalescence in g gens)= (1- 1/N)^g which converges to e^-t
in coalescece units, not time
start with candidate species tree and evaluate
concatination assumes all genes follow same tree-like history

4/27 other topics
distance an d parsimony made be tough on rly divergent sequences
downloaded mrbayes and beast
need to create nexus file for beast and iqtree


