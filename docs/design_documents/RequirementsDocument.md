# Requirements Document

- [Requirements Document](#requirements-document)
- [Overview](#overview)
    - [Figure 1](#figure-1)
- [1. User Requirements](#1-user-requirements)
  - [1.1 User Characteristics](#11-user-characteristics)
  - [1.3 User Interfaces](#13-user-interfaces)
  - [1.4 Functional Requirements](#14-functional-requirements)
    - [**USR.FR.1** Calculate the `AUC` metric for each patient-inhibitor](#usrfr1-calculate-the-auc-metric-for-each-patient-inhibitor)
    - [USR.FR.2 Provide normalization and data cleaning as described by `beatAML-document-placeholder`](#usrfr2-provide-normalization-and-data-cleaning-as-described-by-beataml-document-placeholder)
    - [USR.FR.3 Provide output data in non-technically accessible format](#usrfr3-provide-output-data-in-non-technically-accessible-format)
  - [Non-Functional Requirements](#non-functional-requirements)
    - [USR.NFR.1 Running software should require minimal technical skills](#usrnfr1-running-software-should-require-minimal-technical-skills)
    - [USR.NFR.2 Software errors should prevent data output](#usrnfr2-software-errors-should-prevent-data-output)
    - [USR.NFR.3 ---](#usrnfr3)
- [2. System Requirements](#2-system-requirements)
  - [2.1 Functional Requirements](#21-functional-requirements)
    - [SYS.FR.1 This software must be computationally efficient and parrellizable](#sysfr1-this-software-must-be-computationally-efficient-and-parrellizable)
    - [SYS.FR.2 This software must be well annotated and accessible for maintinience and interative development](#sysfr2-this-software-must-be-well-annotated-and-accessible-for-maintinience-and-interative-development)
  - [2.2 Non-Functional Requirements](#22-non-functional-requirements)
    - [SYS.NFR.1 This software shall `correctly` map (patient, inhbitor) to functional response](#sysnfr1-this-software-shall-correctly-map-patient-inhbitor-to-functional-response)
    - [SYS.NFR.2 This software shall `correctly` calculate AUC metrics.](#sysnfr2-this-software-shall-correctly-calculate-auc-metrics)
    - [SYS.NFR.3 Normalization methods should be `correctly` applied](#sysnfr3-normalization-methods-should-be-correctly-applied)
    - [SYS.NFR.4 Within- & Across- Plate replicates should be handled correctly.](#sysnfr4-within---across--plate-replicates-should-be-handled-correctly)
    - [SYS.NFR.1 This software will be implemented in python and R.](#sysnfr1-this-software-will-be-implemented-in-python-and-r)

# Overview 

This project outlines the software required to process colorimetric `drug response` (MTS assays) data. Data generation follows the procedure: 
1. Patient tumor samples (HNSCC) are obtained from surgery and frozen prior to preparation.  
2. Samples are unfrozen and cultured in the absence of drug 
3. Cells are plated (96-well plate) and drugs are applied. 
4. Cell growth is allowed for ~48 hours 
5. MTS or MTT colorimetric dyes are added (metabolic dye)
6. Optical density is measured (higher values = more cells dead)
7. Results are stored as an excel file in a matrix format corresponding to 96-well plate location, where values correspond to optical density. 

The data processing pipeline that this project comprises requires taking the output of the procedure above and mapping the data to long data format that links patient id, inhibitor, and measured optical density. The full data pipeline is outlined in figure 1: 

### Figure 1

<p align="center">
<img  src="../../figs/pipeline_overview.PNG" >
</p>

Each inhibitor-patient assay has functional response measured over 7 concentrations of inhibitor. 

The output of our data pipeline is the `AUC` metric, which corresponds to patient tumor response to a drug, and is defined as the area under the dose-response curve. 

The ultimate goal of this document is to outline the requirements necessary to develop a robust testing plan to verify and validate the above data pipeline. 

---

# 1. User Requirements

## 1.1 User Characteristics  

There are two user groups: 

1. Lab techs / wet-lab team 
   -  These users have limited programming experience 
   -  Interact only with inputs and outputs of the data but may be involved in testing, identifying bugs, and reading the log files for issues. 
  
2. Bioinformaticians 
     - Highly technical but not necessarily experienced with this package/project. 
     - Most likely will be involved in an advisory capacity with limited time and effort for code review or bug management 
      
## 1.3 User Interfaces   

Users will interact with the data in two formats. The results of the data will be saved to a csv file. This will likely be used in various downstream analysis. Additionally, to facilitate easy visualization, interogation and access there will be a graphical user interface that displays various plots and filtered downloads of the data.

## 1.4 Functional Requirements 

This pipeline was already developed prior to writing this requirements document. 


### **USR.FR.1** Calculate the `AUC` metric for each patient-inhibitor

### USR.FR.2 Provide normalization and data cleaning as described by `beatAML-document-placeholder`

### USR.FR.3 Provide output data in non-technically accessible format

## Non-Functional Requirements

### USR.NFR.1 Running software should require minimal technical skills 

### USR.NFR.2 Software errors should prevent data output 

To avoid non-technical users using faulty output data. 

### USR.NFR.3 --- 

---

# 2. System Requirements  

## 2.1 Functional Requirements

### SYS.FR.1 This software must be computationally efficient and parrellizable 

Assays should be able to process in reasonable time, no slower than 10 seconds per assay. Parrellization, while not necessarily necessary must be optional in order to handle data scaling. 

### SYS.FR.2 This software must be well annotated and accessible for maintinience and interative development 

## 2.2 Non-Functional Requirements

### SYS.NFR.1 This software shall `correctly` map (patient, inhbitor) to functional response

The software will accurately convert photospectrometer matrix-format output to long-data format (lab_id, inhibitor, optical_density) using the specified `plate_maps`.  

### SYS.NFR.2 This software shall `correctly` calculate AUC metrics. 

Area under the dose response curve (AUC) is calculated by integrating a fitted regression (first order, probit regression) the concentration range. AUC values should be deterministic (within minor tolerances) and correctly calculated according to the provided data. 

### SYS.NFR.3 Normalization methods should be `correctly` applied 

These are outlined in `place-holder-normalization` 

### SYS.NFR.4 Within- & Across- Plate replicates should be handled correctly. 

These are outlined in `place-holder-processing` 

### SYS.NFR.1 This software will be implemented in python and R. 

Python for the data pipeline. R shiny for data visualization and interaction in a GUI. 

