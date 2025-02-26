# Rifampicin-PopPK-Repository
## A repository of Rifampicin PopPK models based on MrgSolve by R.
### Overview
This repository contains a collection of Population Pharmacokinetic (PopPK) models for Rifampin, implemented in R. The models are based on various published studies and are designed to simulate pharmacokinetic profiles for both adults and children. The repository includes scripts for model simulation, parameter estimation, and visualization of results. We also offer a easy use of shiny app for our repository, and welcome to download and give us more suggestions.

### Model List
The repository includes 30 different PopPK models for Rifampin, each based on a specific study. Each model is implemented using the mrgsolve and rxode2 packages in R, and the simulations are performed using the PopED package for population pharmacokinetic modeling.

### Repository Structure
simulation_RIF1017.R: The main R script containing the implementation of all 30 Rifampin PopPK models.
README.md: This file, providing an overview of the repository.

### Dependencies
To run the scripts in this repository, you will need the following R packages:
tidyverse: For data manipulation and visualization.
mrgsolve: For pharmacokinetic model simulation.
rxode2: For complex pharmacokinetic model simulation.
cowplot: For combining plots.
ggplot2: For creating visualizations.
ggpubr: For enhanced ggplot2 functionality.
dplyr: For data manipulation.
PopED: For population pharmacokinetic modeling.

### Contributing
Contributions to this repository are welcome! If you would like to add new models, improve existing ones, or fix issues, please follow these steps:
i. Fork the repository.
ii. Create a new branch for your changes: git checkout -b feature/your-feature-name
iii. Commit your changes: git commit -m "Add your commit message here"
iv. Push your changes to your fork: git push origin feature/your-feature-name
v. Open a Pull Request to merge your changes into the main repository.

### License
This project is licensed under the MIT License. See the LICENSE file for details.

## CONTACT
For questions or suggestions, please contact:
Author: Gehang Ju / Graham Ju
Email: 218101055@csu.edu.cn
Institution: Central South University, Xiangya Hospital, Hunan Province, China

## Acknowledgments
This repository is based on published pharmacokinetic studies. The authors of the original studies are acknowledged for their contributions to the field of pharmacokinetics.
