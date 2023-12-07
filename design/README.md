| Tool         	| Purpose                                                   	| Inputs                              	| Outputs                                                     	| Additional parameters                                   	| Scripts? |
|--------------	|-----------------------------------------------------------	|-------------------------------------	|-------------------------------------------------------------	|---------------------------------------------------------	|----------|
| Rosetta      	| "Physics-based" general modeling suite                    	| Structure                           	| Designed sequence, structure with rotamers, energies, stats 	| constraints, restraints, params files, resfile, ...     	|          |
| AlphaFold    	| ML structure prediction                                   	| Sequence                            	| Predicted structure, confidence metrics                     	| templates, MSA, # recycles, # seeds, # models, ...      	| ➕       |
| Protein MPNN 	| ML sequence design                                        	| Protein Structure                   	| Designed sequence                                           	| fixed positions, temperature, ...                       	|          |
| COMBS        	| Statistics-based binding site design                      	| Ligand conformer and template       	| Designed binding site, energies, stats                      	| ligand.txt file, params file, resfile, weightings, ...  	|          |
| RF Diffusion 	| ML structure generation                                   	| Target motif                        	| Generated backbone-only structure                           	| sequence mask, # steps, noise scale, model weights, ... 	| ✅       |
| Reduce       	| Add hydrogens                                             	| Unprotonated structure              	| Protonated structure                                        	|                                                         	| ✅       |
| Probe        	| Find contacts between atoms: vdW, hbonds, clashes         	| All-atom structure                  	| Contacts dataframe                                          	|                                                         	| ✅       |
| "Bunsalyze"  	| Determine buried unsatisfied atoms and clashes in a model 	| All-atom protein + ligand structure 	| Clashing residues, buns atoms, buns energies, etc.          	| Which atoms need hbonds, buns energy weighting          	| ✅       |

This directory includes python functions to run each of these tools, including
* `prepare_xxx.py`: preparing any necessary input (e.g. fixed design positions for protein MPNN)
* `run_xxx.py`: running the tool
* `analyze_xxx.py`: extracting useful statistics from the output (e.g. calculating rmsd of alphafolded structure)

It also includes examples on ways to parallelize these scripts. See our guide to [parallelizing scripts](https://github.com/polizzilab/software-wiki/wiki/Guide:-parallelizing-scripts) for more details:
* with python's `multiprocessing.Pool`
* with bash's `xargs -p`
* and with SLURM

Finally, assorted utility functions can be found under `utils/`.
