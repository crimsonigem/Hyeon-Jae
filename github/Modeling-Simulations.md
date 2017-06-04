## Modeling and Simulations

##### Small Systems Analysis
The production pathways of the different biomaterials that we will be testing (Curli Fibers, Spider Silk, Elastin, Collagen, resilin, Cellulose, Chitosan, Alginate) can be modeled using the laws of Mass Action and Michaelis-Menten Kinetics (for enzyme catalyzed reactions).

This will involve identifying all reactions and interactions involved in the pathway, as well as finding the corresponding reaction coefficients in literature. Once we have our system of differential equations, these can be solved using any computing program.

The 2015 Oxford iGEM team provides a good tutorial on this topic: http://2015.igem.org/Team:Oxford/Modeling/Tutorial

##### Large Systems Modeling (more ambitious)
Genome scale modeling can be used to understand cell responses to metabolic burden.

The iJO1366 model of the K-12 MG1655 *E. coli* strain created by [Orth et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/) can be used as a starting point for generating a genome-scale reconstruction of the strain relevant to our project, LSR10.

Relevant tools in conducting a flux balance analysis are Cobrapy and Escher. Furthermore, the pathways involved in our project can be visualized by importing the genome scale (COBRA) model into EXPASY.

#### References

###### Reaction Kinetics
Grima, Ramon, and Santiago Schnell. [“Modelling Reaction Kinetics inside Cells.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2737326/) Essays in Biochemistry, vol. 45, 2008, pp. 41–56.

###### Single Cell Models  
Cheong, Raymond, et al. [“Models at the Single Cell Level.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3895449/) Wiley Interdisciplinary Reviews: Systems Biology and Medicine, vol. 2, no. 1, 2010, pp. 34–48.

###### Using Modeling Tools  
Jan Schellenberger, et al. [“Quantitative Prediction of Cellular Metabolism with Constraint-Based Models: the COBRA Toolbox v2.0.”](http://www.nature.com/nprot/journal/v6/n9/full/nprot.2011.308.html) Nature Protocols, vol. 6, no. 9, 2011, pp. 1290–307.

###### Genome Scale Modeling  
Archer, Colin T, et al. [“The Genome Sequence of E. Coli W (ATCC 9637): Comparative Genome Analysis and an Improved Genome-Scale Reconstruction of E. Coli.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3032704/) BMC Genomics, vol. 12, 2011, p. 9.

Ines Thiele, and Bernhard Ø Palsson. [“A Protocol for Generating a High-Quality Genome-Scale Metabolic Reconstruction.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125167/) Nature Protocols, vol. 5, no. 1, 2010, pp. 93–121.

O’brien, et al. [“Using Genome-Scale Models to Predict Biological Capabilities.”](https://www-ncbi-nlm-nih-gov/pmc/articles/PMC4451052/) Cell, vol. 161, no. 5, 2015, pp. 971–987.

Orth, Jeffrey D, et al. [“A Comprehensive Genome‐Scale Reconstruction of Escherichia Coli Metabolism—2011.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/) Molecular Systems Biology, vol. 7, no. 1, 2011, p. n/a.

###### Other  
Nash, Trevor R., et al. [“Engineering a Functionalized Biofilm-Based Material for Modulating Escherichia Coli’s Effects in the Mammalian Gastrointestinal Tract.”](https://dash.harvard.edu/bitstream/handle/1/17417585/NASH-SENIORTHESIS-2015.pdf?sequence=1) Engineering a Functionalized Biofilm-Based Material for Modulating Escherichia Coli’s Effects in the Mammalian Gastrointestinal Tract., 2015.

Barnhart, Michelle M., and Chapman, Matthew R. [“Curli Biogenesis and Function.”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2838481/) Annual Review of Microbiology, vol. 60, 2006, p. 131.

###### Past iGEM Projects  
[Imperial College 2016](http://2016.igem.org/Team:Imperial_College/SingleCell) (Co-Cultures)  
[Oxford 2015](http://2015.igem.org/Team:Oxford/Modeling) (UTI Therapeutic w/ AMPs)
