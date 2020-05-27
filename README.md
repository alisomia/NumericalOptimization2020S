## Numerical Optimization Lab 3: Nonnegative Matrix Factorization

Author: Ting Lin @ PKU (lintingsms@pku.edu.cn)

#### Introduction

See `report/main.pdf`

For configuration setting, please run `run_me_first.m` in MATLAB.

#### Folders and Files

```html
./
./README.md							- README file
./run_me_first.m					- Set configuration
./test_syn.n						- test for synthesis data (edit to test)
./nmf_lm_admm2.m					- modified LMF-ADMM
|database/   						- store the database
	V1/								- Synthesis data: case 1
	V2/								- Synthesis data: case 2
	V3/								- Synthesis data: case 3
	V4/								- Synthesis data: case 4
	V5/								- Synthesis data: case 5
	orl.mat							- ORL database (used)
	orl64.mat						- ORL64 database
	yale.mat						- YALE database
	yale64.mat						- YALE64 database (used)
|report/						    - report folder, see `main.pdf`
|src/								- Source code
	|subprob/						- subproblem solver
		|anls/						- solver for ANLS subproblem
			nnls1_asgivens.m		- group activeset solver for ANLS
			nnlsm_activeset.m		- Active Set solver for ANLS
			nnlsm_blockpivot.m		- block pivoting solver for ANLS
			normalEqComb.m			
		|ao_admm/					- solver for AO-ADMM subproblem
			admm_kl_update.m		- KL version
			admm_ls_update.m		- Euclidean version
			terminate.m				- terminate rule
		|apbb/						- solver for APBB subproblem
			pbbnls.m				- projected BB method
		|lm_admm/					- solver for LMF-ADMM subproblem
			JtJprho.m				- LM equation solver
			lm_admm_update.m		- ADMM solver 
			lm_admm_update2.m		
			vec2WH.m				
			WH2vec.m
		|pgd/						- solver for ALSPGD
			nlssubprob.m			- projected gradient descent
		|utils/						- some tools 
			draw.m/					- draw figure by info(s)
			metric.m/				- compute the error
			metric_euc.m/			
			metric_kl.m/
		nmf_admm_euc.m
		nmf_admm_kl.m
		nmf_als.m
		nmf_alspgd.m
		nmf_anls.m
		nmf_anls_activeset.m
		nmf_anls_asgivens.m
		nmf_anls_blockpivot.m
		nmf_ao_admm_euc.m
		nmf_ao_admm_kl.m
		nmf_apbb.m
		nmf_hals.m
		nmf_halsacc.m
		nmf_lm_admm.m
		nmf_mu.m
		nmf_mumod.m
		nmf_muacc.m
		mnf_pgd.m
```





