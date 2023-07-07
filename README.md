# ecFactory: A multi-step method for prediction of metabolic engineering gene targets

![fig1](https://user-images.githubusercontent.com/26483972/175298382-c3ceb172-bf59-4eb7-85e6-aa0ec60928c9.jpg)

The ecFactory method is a series of sequential steps for identification of metabolic engineering gene targets. These targets show which genes should be subject to overexpression, modulated expression (knock-down) or deletion (knock-out), with the objective of increasing production of a given metabolite. This method was developed by combining the principles of the FSEOF algorithm (flux scanning with enforced objective function) together with the features of GECKO enzyme-constrained metabolic models (ecModels), which incorporate enzymes as part of genome-scale metabolic networks.

## Citation
    - 

## Getting started

### Required software
* A functional Matlab installation (MATLAB 7.3 or higher). 
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN.git).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.
* [GECKO 3.0 or higher](https://github.com/SysBioChalmers/GECKO.git)

### Installation

Clone this repository into an accessible directory in your computer. No further steps are needed.

## Tutorial

A case study for prediction of metabolic engineering targets for increased production of 2-phenylethanol in *S. cerevisiae* cells using [ecYeastGEM](https://github.com/SysBioChalmers/yeast-GEM) and the **ecFactory** method is explained in detail in a MATLAB live script. To run this example, open the [live script](https://github.com/SysBioChalmers/ecFactory-case-studies/blob/main/code/find_gene_targets.mlx) in MATLAB and run it! with this, you will see the outputs of the method scripts in real time. 

* An additional case study for prediction of gene targets for enhanced heme production in *S. cerevisiae* has been added. Validation of a subset of the predicted gene targets can be seen [in this publication](https://doi.org/10.1073/pnas.2108245119).

All the relevant outputs of the method are stored in the `tutorials/results` folder in this repository.

## Previous version

Due to significant refactoring of the code, GECKO version 3 is largely not backwards compatible with earlier GECKO versions. Then, compatibility of this repo with initial version of ecFactory is not supported

