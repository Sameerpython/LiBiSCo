# Bindingdata
Ligand Binding Site Comparison (LiBiSCo) is a web-based protein-ligand interaction application for comparing the amino acid residues interacting with atoms of a ligand molecule between different protein-ligand complexes available in the Protein Data Bank (PDB) database. The comparison is performed at the ligand atom level irrespectively of having binding site similarity or not between the protein structures of interest.

The input used in LiBiSCo is a PDB ID or several IDs of a protein-ligand complex(es) and the tool returns a list of identified interactions at ligand atom level including both bonded and non-bonded interactions. A sequence profile for the interaction for each ligand atoms is provided as a WebLogo. The LiBiSco is useful in understanding ligand binding specificity and structural promiscuity among families that are structurally unrelated. 

How to run the server:

1)Install XAMPP

2)Create a folder (For eg. LiBiSCo) inside htdoc folder within XAMPP.

3)Place all the python scripts and image files within LiBiSCo folder.


File details:

The HomePage.py provides two options for the user to analyse their protein-ligand interactions.
	a) Ligand Based Analysis
        b) Cofactor Substructure Based Analysis

a) Ligand Based Analysis:

First the user need to enter PDB codes (for eg. 3WXB,3P19) in the text box provided in the Home page (HomePage.py).  

Second, the user selects the ligand of interest among the list of bound ligands provided and submit for the subsequent analysis (LigPage.py). 

Finally, the results will be provided as tables for ligand–residue interaction (both bonded and non-bonded), color-coded interaction based on physico-chemical properties, and WebLogo for amino acids bound to each ligand atom (Interaction_Result.py).


b) Cofactor Substructure Based Analysis:

First, from the drop-down menu, the user need to select the ligand (cofactor) of interest (HomePage.py). For eg. Select NAD from the dropdown menu and press submit button.

Second, the user need to enter PDB codes (for eg. 3WXB,3P19) in the text box provided in the NAD.py.

Third, the user selects ligand of interest among the list of bound ligands provided for each structure and submit for the subsequent analysis (LigPage1.py).

Finally, the results will be provided as tables for each of the substructure of the selected ligand. For every substructure, ligand–residue interaction (both bonded and non-bonded), color-coded interaction based on physico-chemical properties, and WebLogo for amino acids bound to each ligand atom will be provided (NAD_Results.py).

# Run in Docker

* docker build -t bindingdata .
* docker run -p 80:80 bindingdata:latest
* docker run -p 80:80 -ti bindingdata:latest
