#!/usr/local/Anaconda2.7/bin/python2.7

# Import modules for CGI handling 
import cgi, cgitb
import webbrowser
import urllib
import urllib2
import re
import sys, os
from itertools import izip
import requests
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib
import matplotlib.pyplot as plt
import uuid
from Bio.Seq import Seq
from Bio import motifs
import string
import xlsxwriter
from Bio.Alphabet import IUPAC
from zipfile import ZipFile

# Create instance of FieldStorage
form = cgi.FieldStorage()


print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"
print "<style>"
print """
* {
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  box-sizing: border-box;
}

.grid {
  background: white;
  margin: 0 0 20px 0;
}
.grid:after {
  /* Or @extend clearfix */
  content: "";
  display: table;
  clear: both;
}

[class*='col-'] {
  float: left;
  padding-right: 20px;
}
.grid [class*='col-']:last-of-type {
  padding-right: 0;
}

.col-2-3 {
  width: 33.33%;
  overflow: scroll;
}

.col-1-3 {
  width: 33.33%;
  overflow: scroll;
}



.module {
  padding: 20px;
  background: #eee;
}


body {
  padding: 10px 50px 200px;
  background-size: 300px 300px;
}

h1 {
  color: white;
}
h1 em {
  color: #666;
  font-size: 16px;
}

    """

#############
#style for printing image side by side
print """


.weblogo_column {
    float: left;
    width: 33.33%;
    padding: 2px;
}

/* Clearfix (clear floats) */
.weblogo_row::after {
    content: "";
    clear: both;
    display: table;
}
"""   
############
######style for collapsible content###
print """
.collapsible {
    background-color: #777;
    color: white;
    cursor: pointer;
    padding: 18px;
    width: 100%;
    border: none;
    text-align: left;
    outline: none;
    font-size: 15px;
}

.active, .collapsible:hover {
    background-color: #555;
}

.contentsection {
    padding: 0 18px;
    display: none;
    overflow: scroll;
    background-color: #f1f1f1;
}


"""

############End of style for collapsible content###
#style for divindg into 2 columns

print "* {box-sizing: border-box;}"
print ".column {float: left;width: 50%;padding: 10px;height: 300px;}"
print ".row:after {content: "";display: table;clear: both;}"
#END of style for divindg into 2 columns

print "</style>"
print "<body>"
# Substructure Atom Information for PCHILIDE ligand

Nicot=sorted(['O2','C2','N1','C10','C4X','C4','N3','O4','N10','C9A','C5X','N5', 'C9','C8','C7','C6','C8M','C7M'])
Ribitol=sorted(["O5'","C5'","C4'","O4'","C3'","O3'","C2'", "02'", "C1'"])
phosphate=sorted(['PA','O2A','O1A','O3P','P','O1P','O2P'])
Ribose=sorted(['O5B','C5B','C4B','O4B','C1B','C2B','O2B','C3B','O3B'])

Adenin=sorted(['N6','N1','C2','N3','C4','C5','C6','N7','C8','N9'])

#SUbstructure section ends here

# Information of the selected ligands and PDB ids from LigPage.py
variable = ""
value = ""
r = ""
value_dict={}
lig_sel=[]
for key in form.keys():
        variable = str(key)
#        print "The selected Ligand for PDBID:%s" %variable
        value = str(form.getvalue(variable))
#        print "is", value
        value_dict.setdefault('%s'%variable,[]).append(value)
        r += "<p>"+ variable +", "+ value +"</p>\n"
print "<p style='font-size:20px; color:blue'>  Results for the following selected PDBID's and Ligands: ",'\n'.join("{}:{}".format(k,v) for k,v in value_dict.items()),"</p>","<br/>"


pdbsum_URL="http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode="
pdbsum_URL2="&template=links.html"

#DIctionary and List
pdbsum_dict={}
PDBID_LIST=[]

#Title for Page

Title="The Results are for the following selected PDBID's and Ligands:"
#print Title.center(100,' '),"<br/>"

#Preparing PdbSum Url with selected PDB ids
for ids,lig in value_dict.iteritems():
        pdbsumurl=pdbsum_URL+ids+pdbsum_URL2
        lig_sel.append(lig)
        pdbsum_dict.setdefault('%s'%ids,[]).append(pdbsumurl)
#print "PDBSUMDICT", pdbsum_dict,"<br/>","<br/>"

#creating a list for PDB ids
for id,url in pdbsum_dict.iteritems():
        PDBID_LIST.append(id)
#print "PDBID LIST", PDBID_LIST, "<br/>"

#Extracting the Href links from PDBSum home page for the selected PDB ids and Ligands using BeautifulSoup

litems=[]
new=[]
items2=[]
new1=[]
lig_link=[]
finalLIG_link=[]
liginte_set=set()
ligintelist=[]
link_set=set()

for id,url in pdbsum_dict.iteritems():
    #print id, url, "<br/>"
    for link in url:
        html_page=requests.get(link)
        soup = BeautifulSoup(html_page.text,'html.parser')
        ligand_name_items=soup.find_all('a')
        for items in ligand_name_items:
            name=items.contents[0]
            links='www.ebi.ac.uk' + items.get('href')
            text=str(name)+ " " + str(links)
            litems.append(text)

#Looping over the extracted URLs from PDBSum and appending into a lIst called New
for x in litems:
        x=x.strip()
        new.append(x)
#print new

#Looping over Ligand and PDBSum Urls (from above step) to extract the PDBSUm URL for the seleted Ligand Page in PDBSUm. The Ligand Page URL is now as a SET data type
ligand_urlLIST=[]
for lig in lig_sel:
    #print "LIGAND", lig, "<br/>"
    for y in new:
        lig= ''.join(lig)
        if y.startswith(lig):#y=y.split()
            y=y.split()
            link=y[1]
            link1="http://"+link
            #print link1
            ligand_urlLIST.append(link1)
            link_set.add(link1)
    link_setlist=list(link_set)
#print ligand_urlLIST

#print "SET",link_setlist, "<br/>"

#Using Beautifulsoup to extract all the links from PDBSUM Ligand interaction Page for each of the PDB ids.
#print PDBID_LIST
PDBID_URL_dict=zip(PDBID_LIST,ligand_urlLIST)
LiginteractPage=dict(PDBID_URL_dict)
#print "what1", LiginteractPage, "<br/>"
#print "what2", PDBID_URL_dict

for pdbid,pdbsumlink in LiginteractPage.iteritems():
        links=list(LiginteractPage.viewvalues())
        for y in links:
                html_page2=requests.get(y)
                soup2 = BeautifulSoup(html_page2.text,'html.parser')
                ligand_name_items1=soup2.find_all('a')

#Looping over all the href links from PDBSUM ligand intercation page to extract 
#URL(Final page for extracting atom details) for the atom based interaction for PDB ids with ligands

                for items in ligand_name_items1:
                        name=items.contents[0]
                        links='www.ebi.ac.uk' + items.get('href')
                        text=str(links)
                        items2.append(text)

                for i in items2:
#                       print i
                        final=i.split()
                        final=''.join(final)
                        #print final, "<br/>"
                        finalLIG_link.append(final)
                        lastitem=finalLIG_link[-1]
                        lastitem="http://"+lastitem
#                print lastitem,"<br/>"

                if lastitem not in ligintelist:
                        ligintelist.append(lastitem)
                liginte_set.add(lastitem)
        liginte_list=list(liginte_set)
#print "LIST", (ligintelist),"<br/>"
PDBID_INTURL_dict=zip(PDBID_LIST,ligintelist)
pdbsum_dict1=dict(PDBID_INTURL_dict)
for id,link in pdbsum_dict1.iteritems():
        links=list(pdbsum_dict1.viewvalues())
        PDBID=list(pdbsum_dict1.viewkeys())
#print "HI LINKS CHECK for SPACE", links,"<br/>"
#print "PDBIDs", PDBID,"<br/>"

Number_of_Ids=len(PDBID)
#FInal DIctionary with PDBID and Ligplot URL for extracting intercation details
mydictcheck={}
for ids,links in zip(PDBID,ligintelist):
        mydictcheck.setdefault('%s'%ids,[]).append(links)

###############################################################################
        #SECTION OF FIDING COMMON LIGAND ATOMS

###############################################################################       

#selecting common ligand atoms that are hydrogen bonded in selected PDB structures
H_printing = False
H_atoms_commoncomp={}
for H_pdbids,H_pdbidlinks in mydictcheck.iteritems():
    for H_links_sel in H_pdbidlinks:
        H_links_sel1=str(H_links_sel)
        weblink=requests.get(H_links_sel1, stream=True)
        for H_atomlines in weblink.iter_lines():
            H_atomlines1=H_atomlines.strip()
            
            if H_atomlines1.startswith('Hydrogen bonds'):
                H_printing = True
            elif H_atomlines1.startswith('Non-bonded contacts'):
                H_printing = False
            if H_printing:
                #print H_atomlines
                if H_atomlines1.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                    H_atomlines2=H_atomlines1.split()
                    H_atm_sel=H_atomlines2[8]
                    H_atoms_commoncomp.setdefault('%s'%H_pdbids,[]).append(H_atm_sel)

H_atomsvalues_dict1=H_atoms_commoncomp.values()
H_common_intersectionfinal=sorted(list(set.intersection(*map(set,H_atomsvalues_dict1))))

#END of selecting common ligand atoms that are hydrogen bonded in selected PDB structures

#selecting common ligand atoms that are non-hydrogen bonded in selected PDB structures
NONH_printing = False
NONHatoms_commoncomp={}
for NONHpdbids,NONHpdbidlinks in mydictcheck.iteritems():
    
    for NONHlinks_sel in NONHpdbidlinks:
        NONHlinks_sel1=str(NONHlinks_sel)
        weblink=requests.get(NONHlinks_sel1, stream=True)
        for NONHatomlines in weblink.iter_lines():
            NONHatomlines1=NONHatomlines.strip()
            if NONHatomlines1.startswith('Non-bonded contacts'):
                NONH_printing = True
                
            elif NONHatomlines1.startswith('Hydrogen bonds'):
                NONH_printing = False
            if NONH_printing:                
                #print atomlines1
                if NONHatomlines1.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                    
       
                    NONHatomlines2=NONHatomlines1.split()
                    NONHatm_sel=NONHatomlines2[8]
                    NONHatoms_commoncomp.setdefault('%s'%NONHpdbids,[]).append(NONHatm_sel)

NONHatomsvalues_dict1=NONHatoms_commoncomp.values()
NONHcommon_intersectionfinal=sorted(list(set.intersection(*map(set,NONHatomsvalues_dict1))))                

###############################################################################
        #END OF SECTION OF FIDING COMMON LIGAND ATOMS

############################################################################### 

###############################################################################
        #START OF SECTION OF LIGAND ATOMS IN EACH SUBGROUPS:

############################################################################### 

###############################################################################
        # 1.START OF SECTION OF NAD SUBGROUPS:

############################################################################### 
#Nicot
Nicot_graphdicH={}
Nicot_common_graphdicH={}
Nicot_graphdicNH={}
Nicot_common_graphdicNH={}

Nicot_All_combine_Lig_Res_H={}
Nicot_allH_Lig_Resdict={}
Nicot_Common_combine_Lig_Res_H={}
Nicot_CommonH_Lig_Resdict={}

Nicot_All_combine_Lig_Res_NH={}
Nicot_allNH_Lig_Resdict={}
Nicot_Common_combine_Lig_Res_NH={}
Nicot_CommonNH_Lig_Resdict={}

Nicot_All_combine_Lig_Res_H_distance={}
Nicot_allH_Lig_Resdict_distance={}
Nicot_Common_combine_Lig_Res_H_distance={}
Nicot_CommonH_Lig_Resdict_distance={}


Nicot_All_combine_Lig_Res_NH_distance={}
Nicot_allNH_Lig_Resdict_distance={}
Nicot_Common_combine_Lig_Res_NH_distance={}
Nicot_CommonNH_Lig_Resdict_distance={}


Nicot_listdata_H=[]
Nicot_listdata_NH=[]

Nicot_lresidueH=[]
Nicot_latomH=[]
Nicot_lresidueNH=[]
Nicot_latomNH=[]

Nicot_common_listdata_H=[]
Nicot_common_listdata_NH=[]
Nicot_H_appended_lig_tabledic={}
Nicot_H_Common_appended_lig_tabledic={}

Nicot_NH_appended_lig_tabledic={}
Nicot_NH_Common_appended_lig_tabledic={}
Nicot_finalsetH=set()
Nicot_finalsetNH=set()
#End of Nicot
#Ribitol
Ribitol_graphdicH={}
Ribitol_common_graphdicH={}
Ribitol_graphdicNH={}
Ribitol_common_graphdicNH={}

Ribitol_All_combine_Lig_Res_H={}
Ribitol_allH_Lig_Resdict={}
Ribitol_Common_combine_Lig_Res_H={}
Ribitol_CommonH_Lig_Resdict={}

Ribitol_All_combine_Lig_Res_NH={}
Ribitol_allNH_Lig_Resdict={}
Ribitol_Common_combine_Lig_Res_NH={}
Ribitol_CommonNH_Lig_Resdict={}

Ribitol_All_combine_Lig_Res_H_distance={}
Ribitol_allH_Lig_Resdict_distance={}
Ribitol_Common_combine_Lig_Res_H_distance={}
Ribitol_CommonH_Lig_Resdict_distance={}


Ribitol_All_combine_Lig_Res_NH_distance={}
Ribitol_allNH_Lig_Resdict_distance={}
Ribitol_Common_combine_Lig_Res_NH_distance={}
Ribitol_CommonNH_Lig_Resdict_distance={}



Ribitol_listdata_H=[]
Ribitol_listdata_NH=[]

Ribitol_lresidueH=[]
Ribitol_latomH=[]
Ribitol_lresidueNH=[]
Ribitol_latomNH=[]

Ribitol_common_listdata_H=[]
Ribitol_common_listdata_NH=[]
Ribitol_H_appended_lig_tabledic={}
Ribitol_H_Common_appended_lig_tabledic={}

Ribitol_NH_appended_lig_tabledic={}
Ribitol_NH_Common_appended_lig_tabledic={}
#End of Ribitol

#phosphate
phosphate_graphdicH={}
phosphate_common_graphdicH={}
phosphate_graphdicNH={}
phosphate_common_graphdicNH={}

phosphate_All_combine_Lig_Res_H={}
phosphate_allH_Lig_Resdict={}
phosphate_Common_combine_Lig_Res_H={}
phosphate_CommonH_Lig_Resdict={}

phosphate_All_combine_Lig_Res_NH={}
phosphate_allNH_Lig_Resdict={}
phosphate_Common_combine_Lig_Res_NH={}
phosphate_CommonNH_Lig_Resdict={}

phosphate_All_combine_Lig_Res_H_distance={}
phosphate_allH_Lig_Resdict_distance={}
phosphate_Common_combine_Lig_Res_H_distance={}
phosphate_CommonH_Lig_Resdict_distance={}


phosphate_All_combine_Lig_Res_NH_distance={}
phosphate_allNH_Lig_Resdict_distance={}
phosphate_Common_combine_Lig_Res_NH_distance={}
phosphate_CommonNH_Lig_Resdict_distance={}


phosphate_listdata_H=[]
phosphate_listdata_NH=[]

phosphate_lresidueH=[]
phosphate_latomH=[]
phosphate_lresidueNH=[]
phosphate_latomNH=[]

phosphate_common_listdata_H=[]
phosphate_common_listdata_NH=[]
phosphate_H_appended_lig_tabledic={}
phosphate_H_Common_appended_lig_tabledic={}

phosphate_NH_appended_lig_tabledic={}
phosphate_NH_Common_appended_lig_tabledic={}
phosphate_finalsetH=set()
phosphate_finalsetNH=set()
#End of phosphate


#Ribose
Ribose_graphdicH={}
Ribose_common_graphdicH={}
Ribose_graphdicNH={}
Ribose_common_graphdicNH={}

Ribose_All_combine_Lig_Res_H={}
Ribose_allH_Lig_Resdict={}
Ribose_Common_combine_Lig_Res_H={}
Ribose_CommonH_Lig_Resdict={}

Ribose_All_combine_Lig_Res_NH={}
Ribose_allNH_Lig_Resdict={}
Ribose_Common_combine_Lig_Res_NH={}
Ribose_CommonNH_Lig_Resdict={}

Ribose_All_combine_Lig_Res_H_distance={}
Ribose_allH_Lig_Resdict_distance={}
Ribose_Common_combine_Lig_Res_H_distance={}
Ribose_CommonH_Lig_Resdict_distance={}


Ribose_All_combine_Lig_Res_NH_distance={}
Ribose_allNH_Lig_Resdict_distance={}
Ribose_Common_combine_Lig_Res_NH_distance={}
Ribose_CommonNH_Lig_Resdict_distance={}


Ribose_listdata_H=[]
Ribose_listdata_NH=[]

Ribose_common_listdata_H=[]
Ribose_common_listdata_NH=[]
Ribose_H_appended_lig_tabledic={}
Ribose_H_Common_appended_lig_tabledic={}

Ribose_NH_appended_lig_tabledic={}
Ribose_NH_Common_appended_lig_tabledic={}
#End ofRibose
#Adenin
Adenin_graphdicH={}
Adenin_common_graphdicH={}
Adenin_graphdicNH={}
Adenin_common_graphdicNH={}

Adenin_All_combine_Lig_Res_H={}
Adenin_allH_Lig_Resdict={}
Adenin_Common_combine_Lig_Res_H={}
Adenin_CommonH_Lig_Resdict={}

Adenin_All_combine_Lig_Res_NH={}
Adenin_allNH_Lig_Resdict={}
Adenin_Common_combine_Lig_Res_NH={}
Adenin_CommonNH_Lig_Resdict={}

Adenin_All_combine_Lig_Res_H_distance={}
Adenin_allH_Lig_Resdict_distance={}
Adenin_Common_combine_Lig_Res_H_distance={}
Adenin_CommonH_Lig_Resdict_distance={}


Adenin_All_combine_Lig_Res_NH_distance={}
Adenin_allNH_Lig_Resdict_distance={}
Adenin_Common_combine_Lig_Res_NH_distance={}
Adenin_CommonNH_Lig_Resdict_distance={}


Adenin_listdata_H=[]
Adenin_listdata_NH=[]

Adenin_lresidueH=[]
Adenin_latomH=[]
Adenin_lresidueNH=[]
Adenin_latomNH=[]

Adenin_common_listdata_H=[]
Adenin_common_listdata_NH=[]
Adenin_H_appended_lig_tabledic={}
Adenin_H_Common_appended_lig_tabledic={}

Adenin_NH_appended_lig_tabledic={}
Adenin_NH_Common_appended_lig_tabledic={}
#End ofAdenin
lresidueH=[]
latomH=[]
ldistanceH=[]
residueH={}
atmnameH={}
dicresidue_unique={}
residue_seenH=set()
atom_seenH=set()

lresidueNH=[]
latomNH=[]
ldistanceNH=[]
residueNH={}
atmnameNH={}
dicresidue_unique={}
residue_seenNH=set()
atom_seenNH=set()

Ribitol_finalsetH=set()
Ribitol_finalsetNH=set()
Adenin_finalsetH=set()
Adenin_finalsetNH=set()

finalsetH=set()
finalsetNH=set()
combines_listdata=[]
#print "<table style=width:50%>"
#print "<tr>"
#print "<th colspan='%d'>Interaction List</th>"% Number_of_Ids

#print "</tr>"
#print "</div>"

H_appended_lig_tabledic={}
H_Common_appended_lig_tabledic={}
NH_Common_appended_lig_tabledic={}
NH_appended_lig_tabledic={}
graphdicH={}
graphdicNH={}

common_graphdicH={}
common_graphdicNH={}
printing = False

for id,link in mydictcheck.iteritems():
   #print link, id	   
   links_sel=link[0]
   link1= ''.join(str(links_sel))
   res2=urllib.urlopen(str(link1))
   html=res2.read()
   #print html
   for l in link:
            ll=str(l)
            r = requests.get(ll, stream=True)
            for line in r.iter_lines():
                    line=line.strip()
                    if line.startswith('Hydrogen bonds'):
                        printing = True
                    elif line.startswith('Non-bonded contacts'):
                        printing = False 
                    if printing:
                        if line.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                            #print "HB", line
                            lineH=line.split()
                            lignameH=lineH[9]
                            atmH=lineH[8]
                            resH=lineH[3]
                            residuenumH=lineH[4]
                            distanceH=lineH[12]
                            resnumH=resH+residuenumH
#                               #appending each residue and its position to list called lresidue
                            lresidueH.append(resnumH)
                            #appending each ligand atom to list called latom
                            latomH.append(atmH)
                            #appending distance of each interaction to ldistance
                            ldistanceH.append(distanceH)
                            #creating a set for residue with position
                            residue_seenH.add(resnumH)
                            #creating a set for each ligand atom
                            atom_seenH.add(atmH)
                            #making a dictionary with list comtaining residue name and position
                            residueH.setdefault('%s'%id,[]).append(resnumH)
                            #making a dictionary with list comataing ligand atoms
                            atmnameH.setdefault('%s'%id,[]).append(atmH)
                            if atmH in Nicot:
                                Nicot_lresidueH.append(resnumH)
                                Nicot_latomH.append(atmH)
                                Nicot_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for physio and weblogo
                                Nicot_All_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for table
                                Nicot_All_combine_Lig_Res_H_uniquify= {k:list(set(j)) for k,j in Nicot_All_combine_Lig_Res_H.items()}
                                Nicot_allH_Lig_Resdict['%s'%id]=Nicot_All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group
                                
                                Nicot_All_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                Nicot_All_combine_Lig_Res_H_distance_uniquify= {k:list(set(j)) for k,j in Nicot_All_combine_Lig_Res_H_distance.items()}
                                Nicot_allH_Lig_Resdict_distance['%s'%id]=Nicot_All_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                                if atmH in H_common_intersectionfinal:
                                    Nicot_common_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)
                                    Nicot_Common_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with common lig atom and residues for table
                                    Nicot_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Nicot_Common_combine_Lig_Res_H.items()}
                                    Nicot_CommonH_Lig_Resdict['%s'%id]=Nicot_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Nicot_Common_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                    Nicot_Common_combine_Lig_Res_H_distance_uniquify={k:list(set(j)) for k,j in Nicot_Common_combine_Lig_Res_H_distance.items()}
                                    Nicot_CommonH_Lig_Resdict_distance['%s'%id]=Nicot_Common_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
 
                            if atmH in Ribitol:
                                Ribitol_lresidueH.append(resnumH)
                                Ribitol_latomH.append(atmH)
                                Ribitol_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for physio and weblogo
                                Ribitol_All_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for table
                                Ribitol_All_combine_Lig_Res_H_uniquify= {k:list(set(j)) for k,j in Ribitol_All_combine_Lig_Res_H.items()}
                                Ribitol_allH_Lig_Resdict['%s'%id]=Ribitol_All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Ribitol_All_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                Ribitol_All_combine_Lig_Res_H_distance_uniquify= {k:list(set(j)) for k,j in Ribitol_All_combine_Lig_Res_H_distance.items()}
                                Ribitol_allH_Lig_Resdict_distance['%s'%id]=Ribitol_All_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                                
                                if atmH in H_common_intersectionfinal:
                                    Ribitol_common_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)
                                    Ribitol_Common_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with common lig atom and residues for table
                                    Ribitol_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Ribitol_Common_combine_Lig_Res_H.items()}
                                    Ribitol_CommonH_Lig_Resdict['%s'%id]=Ribitol_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Ribitol_Common_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                    Ribitol_Common_combine_Lig_Res_H_distance_uniquify={k:list(set(j)) for k,j in Ribitol_Common_combine_Lig_Res_H_distance.items()}
                                    Ribitol_CommonH_Lig_Resdict_distance['%s'%id]=Ribitol_Common_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                            if atmH in phosphate:
                                phosphate_lresidueH.append(resnumH)
                                phosphate_latomH.append(atmH)
                                phosphate_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for physio and weblogo
                                phosphate_All_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for table
                                phosphate_All_combine_Lig_Res_H_uniquify= {k:list(set(j)) for k,j in phosphate_All_combine_Lig_Res_H.items()}
                                phosphate_allH_Lig_Resdict['%s'%id]=phosphate_All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                phosphate_All_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                phosphate_All_combine_Lig_Res_H_distance_uniquify= {k:list(set(j)) for k,j in phosphate_All_combine_Lig_Res_H_distance.items()}
                                phosphate_allH_Lig_Resdict_distance['%s'%id]=phosphate_All_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                                
                                if atmH in H_common_intersectionfinal:
                                    phosphate_common_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)
                                    phosphate_Common_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with common lig atom and residues for table
                                    phosphate_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in phosphate_Common_combine_Lig_Res_H.items()}
                                    phosphate_CommonH_Lig_Resdict['%s'%id]=phosphate_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    phosphate_Common_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                    phosphate_Common_combine_Lig_Res_H_distance_uniquify={k:list(set(j)) for k,j in phosphate_Common_combine_Lig_Res_H_distance.items()}
                                    phosphate_CommonH_Lig_Resdict_distance['%s'%id]=phosphate_Common_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                            if atmH in Ribose:

                                Ribose_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)
                                Ribose_All_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for table
                                Ribose_All_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Ribose_All_combine_Lig_Res_H.items()}
                                Ribose_allH_Lig_Resdict['%s'%id]=Ribose_All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Ribose_All_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                Ribose_All_combine_Lig_Res_H_distance_uniquify= {k:list(set(j)) for k,j in Ribose_All_combine_Lig_Res_H_distance.items()}
                                Ribose_allH_Lig_Resdict_distance['%s'%id]=Ribose_All_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                

                                if atmH in H_common_intersectionfinal:
                                    Ribose_common_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)                                    
                                    Ribose_Common_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with common lig atom and residues for table
                                    Ribose_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Ribose_Common_combine_Lig_Res_H.items()}
                                    Ribose_CommonH_Lig_Resdict['%s'%id]=Ribose_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Ribose_Common_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                    Ribose_Common_combine_Lig_Res_H_distance_uniquify={k:list(set(j)) for k,j in Ribose_Common_combine_Lig_Res_H_distance.items()}
                                    Ribose_CommonH_Lig_Resdict_distance['%s'%id]=Ribose_Common_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                

 
                            if atmH in Adenin:
                                Adenin_lresidueH.append(resnumH)
                                Adenin_latomH.append(atmH)                                
                                Adenin_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)
                                Adenin_All_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with all lig atom and residues for table
                                Adenin_All_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Adenin_All_combine_Lig_Res_H.items()}
                                Adenin_allH_Lig_Resdict['%s'%id]=Adenin_All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Adenin_All_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                Adenin_All_combine_Lig_Res_H_distance_uniquify= {k:list(set(j)) for k,j in Adenin_All_combine_Lig_Res_H_distance.items()}
                                Adenin_allH_Lig_Resdict_distance['%s'%id]=Adenin_All_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                                if atmH in H_common_intersectionfinal:
                                    Adenin_common_graphdicH.setdefault('%s'%atmH,[]).append(resnumH)                                    
                                    Adenin_Common_combine_Lig_Res_H.setdefault('%s'%atmH,[]).append(resnumH)#creating dictionary with common lig atom and residues for table
                                    Adenin_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Adenin_Common_combine_Lig_Res_H.items()}
                                    Adenin_CommonH_Lig_Resdict['%s'%id]=Adenin_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Adenin_Common_combine_Lig_Res_H_distance.setdefault('%s'%atmH,[]).append(distanceH)#creating dictionary with all lig atom and distance for table
                                    Adenin_Common_combine_Lig_Res_H_distance_uniquify={k:list(set(j)) for k,j in Adenin_Common_combine_Lig_Res_H_distance.items()}
                                    Adenin_CommonH_Lig_Resdict_distance['%s'%id]=Adenin_Common_combine_Lig_Res_H_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

                        
                    else:
                        if line.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                            #print "NHB", line
                            lineNH=line.split()
                            lignameNH=lineNH[9]
                            atmNH=lineNH[8]
                            resNH=lineNH[3]
                            residuenumNH=lineNH[4]
                            distanceNH=lineNH[12]
                            resnumNH=resNH+residuenumNH
#                               #appending each residue and its position to list called lresidue
                            lresidueNH.append(resnumNH)
                            #appending each ligand atom to list called latom
                            latomNH.append(atmNH)
                            #appending distance of each interaction to ldistance
                            ldistanceNH.append(distanceNH)
                            #creating a set for residue with position
                            residue_seenNH.add(resnumNH)
                            #creating a set for each ligand atom
                            atom_seenNH.add(atmNH)
                            #making a dictionary with list comtaining residue name and position
                            residueNH.setdefault('%s'%id,[]).append(resnumNH)
                            #making a dictionary with list comataing ligand atoms
                            atmnameNH.setdefault('%s'%id,[]).append(atmNH)
                            if atmNH in Nicot:
                                Nicot_lresidueNH.append(resnumNH)
                                Nicot_latomNH.append(atmNH)
                                Nicot_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)
                                Nicot_All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                Nicot_All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Nicot_All_combine_Lig_Res_NH.items()}
                                Nicot_allNH_Lig_Resdict['%s'%id]=Nicot_All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Nicot_All_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and residues for table
                                Nicot_All_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Nicot_All_combine_Lig_Res_NH_distance.items()}
                                Nicot_allNH_Lig_Resdict_distance['%s'%id]=Nicot_All_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                if atmNH in NONHcommon_intersectionfinal:
                                    Nicot_common_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)                                   
                                    Nicot_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    Nicot_Common_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Nicot_Common_combine_Lig_Res_NH.items()}
                                    Nicot_CommonNH_Lig_Resdict['%s'%id]=Nicot_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Nicot_Common_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with common lig atom and residues for table
                                    Nicot_Common_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Nicot_Common_combine_Lig_Res_NH_distance.items()}
                                    Nicot_CommonNH_Lig_Resdict_distance['%s'%id]=Nicot_Common_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for common group

                            if atmNH in Ribitol:
                                Ribitol_lresidueNH.append(resnumNH)
                                Ribitol_latomNH.append(atmNH)
                                Ribitol_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)
                                Ribitol_All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                Ribitol_All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Ribitol_All_combine_Lig_Res_NH.items()}
                                Ribitol_allNH_Lig_Resdict['%s'%id]=Ribitol_All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Ribitol_All_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and residues for table
                                Ribitol_All_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Ribitol_All_combine_Lig_Res_NH_distance.items()}
                                Ribitol_allNH_Lig_Resdict_distance['%s'%id]=Ribitol_All_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for all group
                                

                                if atmNH in NONHcommon_intersectionfinal:
                                    Ribitol_common_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)                                   
                                    Ribitol_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    Ribitol_Common_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Ribitol_Common_combine_Lig_Res_NH.items()}
                                    Ribitol_CommonNH_Lig_Resdict['%s'%id]=Ribitol_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Ribitol_Common_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with common lig atom and residues for table
                                    Ribitol_Common_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Ribitol_Common_combine_Lig_Res_NH_distance.items()}
                                    Ribitol_CommonNH_Lig_Resdict_distance['%s'%id]=Ribitol_Common_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for common group

                            if atmNH in phosphate:
                                phosphate_lresidueNH.append(resnumNH)
                                phosphate_latomNH.append(atmNH)
                                phosphate_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)
                                phosphate_All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                phosphate_All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in phosphate_All_combine_Lig_Res_NH.items()}
                                phosphate_allNH_Lig_Resdict['%s'%id]=phosphate_All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                phosphate_All_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and residues for table
                                phosphate_All_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in phosphate_All_combine_Lig_Res_NH_distance.items()}
                                phosphate_allNH_Lig_Resdict_distance['%s'%id]=phosphate_All_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for all group                                

                                if atmNH in NONHcommon_intersectionfinal:
                                    phosphate_common_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)                                   
                                    phosphate_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    phosphate_Common_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in phosphate_Common_combine_Lig_Res_NH.items()}
                                    phosphate_CommonNH_Lig_Resdict['%s'%id]=phosphate_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    phosphate_Common_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with common lig atom and residues for table
                                    phosphate_Common_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in phosphate_Common_combine_Lig_Res_NH_distance.items()}
                                    phosphate_CommonNH_Lig_Resdict_distance['%s'%id]=phosphate_Common_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for common group

                            if atmNH in Ribose:
                                Ribose_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)
                                Ribose_All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                Ribose_All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Ribose_All_combine_Lig_Res_NH.items()}
                                Ribose_allNH_Lig_Resdict['%s'%id]=Ribose_All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                Ribose_All_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and residues for table
                                Ribose_All_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Ribose_All_combine_Lig_Res_NH_distance.items()}
                                Ribose_allNH_Lig_Resdict_distance['%s'%id]=Ribose_All_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                if atmNH in NONHcommon_intersectionfinal:
                                    Ribose_common_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)                                   
                                    Ribose_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    Ribose_Common_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Ribose_Common_combine_Lig_Res_NH.items()}
                                    Ribose_CommonNH_Lig_Resdict['%s'%id]=Ribose_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    Ribose_Common_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with common lig atom and residues for table
                                    Ribose_Common_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Ribose_Common_combine_Lig_Res_NH_distance.items()}
                                    Ribose_CommonNH_Lig_Resdict_distance['%s'%id]=Ribose_Common_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for common group

                            if atmNH in Adenin:
                                Adenin_lresidueNH.append(resnumNH)
                                Adenin_latomNH.append(atmNH)
                                Adenin_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)
                                Adenin_All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                Adenin_All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Adenin_All_combine_Lig_Res_NH.items()}
                                Adenin_allNH_Lig_Resdict['%s'%id]=Adenin_All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group
                                
                                


                                Adenin_All_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and residues for table
                                Adenin_All_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Adenin_All_combine_Lig_Res_NH_distance.items()}
                                Adenin_allNH_Lig_Resdict_distance['%s'%id]=Adenin_All_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for all group
                                


                                if atmNH in NONHcommon_intersectionfinal:
                                    Adenin_common_graphdicNH.setdefault('%s'%atmNH,[]).append(resnumNH)                                   
                                    Adenin_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    Adenin_Common_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in Adenin_Common_combine_Lig_Res_NH.items()}
                                    Adenin_CommonNH_Lig_Resdict['%s'%id]=Adenin_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group


                                    Adenin_Common_combine_Lig_Res_NH_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with common lig atom and residues for table
                                    Adenin_Common_combine_Lig_Res_NH_distance_uniquify= {k:list(set(j)) for k,j in Adenin_Common_combine_Lig_Res_NH_distance.items()}
                                    Adenin_CommonNH_Lig_Resdict_distance['%s'%id]=Adenin_Common_combine_Lig_Res_NH_distance_uniquify#final dic for table with pdb id , lig atom and residues for common group


                       
   
   



   Nicot_listdata_H=[]
   Nicot_listdata_NH=[]
   Nicot_lresidueH=[]
   Nicot_latomH=[]
   Nicot_lresidueNH=[]
   Nicot_latomNH=[]
   Nicot_All_combine_Lig_Res_H={}
   Nicot_Common_combine_Lig_Res_H={}
   Nicot_All_combine_Lig_Res_NH={}
   Nicot_Common_combine_Lig_Res_NH={}

   Nicot_All_combine_Lig_Res_H_distance={}
   Nicot_Common_combine_Lig_Res_H_distance={}
   Nicot_All_combine_Lig_Res_NH_distance={}
   Nicot_Common_combine_Lig_Res_NH_distance={}
   
   Nicot_common_listdata_H=[]
   Nicot_common_listdata_NH=[]
   Nicot_finalsetH.clear()
   Nicot_finalsetNH.clear()

   Ribitol_listdata_H=[]
   Ribitol_listdata_NH=[]
   Ribitol_lresidueH=[]
   Ribitol_latomH=[]
   Ribitol_lresidueNH=[]
   Ribitol_latomNH=[]
   Ribitol_All_combine_Lig_Res_H={}
   Ribitol_Common_combine_Lig_Res_H={}
   Ribitol_All_combine_Lig_Res_NH={}
   Ribitol_Common_combine_Lig_Res_NH={}
   Ribitol_All_combine_Lig_Res_H_distance={}
   Ribitol_Common_combine_Lig_Res_H_distance={}
   Ribitol_All_combine_Lig_Res_NH_distance={}
   Ribitol_Common_combine_Lig_Res_NH_distance={}


   phosphate_listdata_H=[]
   phosphate_listdata_NH=[]
   phosphate_lresidueH=[]
   phosphate_latomH=[]
   phosphate_lresidueNH=[]
   phosphate_latomNH=[]
   phosphate_All_combine_Lig_Res_H={}
   phosphate_Common_combine_Lig_Res_H={}
   phosphate_All_combine_Lig_Res_NH={}
   phosphate_Common_combine_Lig_Res_NH={}
   phosphate_All_combine_Lig_Res_H_distance={}
   phosphate_Common_combine_Lig_Res_H_distance={}
   phosphate_All_combine_Lig_Res_NH_distance={}
   phosphate_Common_combine_Lig_Res_NH_distance={}

   phosphate_common_listdata_H=[]
   phosphate_common_listdata_NH=[]
   phosphate_finalsetH.clear()
   phosphate_finalsetNH.clear()


   
   Ribose_listdata_H=[]
   Ribose_listdata_NH=[]
   Ribose_All_combine_Lig_Res_H={}
   Ribose_Common_combine_Lig_Res_H={}   
   Ribose_All_combine_Lig_Res_NH={}
   Ribose_Common_combine_Lig_Res_NH={}   
   Ribose_All_combine_Lig_Res_H_distance={}
   Ribose_Common_combine_Lig_Res_H_distance={}
   Ribose_All_combine_Lig_Res_NH_distance={}
   Ribose_Common_combine_Lig_Res_NH_distance={}

   
   Adenin_lresidueH=[]
   Adenin_latomH=[]
   Adenin_lresidueNH=[]
   Adenin_latomNH=[]
   Adenin_listdata_H=[]
   Adenin_listdata_NH=[]
   Adenin_All_combine_Lig_Res_H={}
   Adenin_Common_combine_Lig_Res_H={}   
   Adenin_All_combine_Lig_Res_NH={}
   Adenin_Common_combine_Lig_Res_NH={}   
   Adenin_All_combine_Lig_Res_H_distance={}
   Adenin_Common_combine_Lig_Res_H_distance={}
   Adenin_All_combine_Lig_Res_NH_distance={}
   Adenin_Common_combine_Lig_Res_NH_distance={}

   
   Ribitol_common_listdata_H=[]
   Ribitol_common_listdata_NH=[]
   Ribose_common_listdata_H=[]
   Ribose_common_listdata_NH=[]
   Adenin_common_listdata_H=[]
   Adenin_common_listdata_NH=[]

   
   lresidueH=[]
   latomH=[]
   ldistanceH=[]
   
   lresidueNH=[]
   latomNH=[]
   ldistanceNH=[]
   
   combines_listdata=[]
   residue_seenH.clear()
   Ribitol_finalsetH.clear()
   Ribitol_finalsetNH.clear()
   Adenin_finalsetH.clear()
   Adenin_finalsetNH.clear()

   finalsetH.clear()
   finalsetNH.clear()

####################Define function for Statistics ################################
def percentage(dictname,subgroup):
    
    Count_Atom={}
    percentage_Atom={}
    atmlist=[]
    if bool(dictname):
        for key, value in dictname.iteritems():
            
            for atom in subgroup:
            
                for key1,value1 in value.iteritems():
                #for i in dict1.keys():
                    if atom == key1:
                        Count_Atom[key1]=1
                        percentage_Atom['%s'%key]=Count_Atom
                        #print percent
            Count_Atom={}
        tabl=pd.DataFrame.from_dict(percentage_Atom).fillna(0)
        Num_cols = len (PDBID_LIST)
        for atms in percentage_Atom.values():
            for atms_key in atms.keys():
                atmlist.append(atms_key)
            count_atmlist=list(set(atmlist))
        tabl['Percentage of Interaction']= (tabl.sum(axis=1)/Num_cols)*100
        tabl['Percentage of Interaction']=tabl['Percentage of Interaction'].round(2)
        print "<br/>"," No. of Ligand atoms:", len(count_atmlist), "/",len(subgroup), "<br/>"
        print tabl.T.to_html(justify='center'),"<br/>"
        #print tabl.style.background_gradient(cmap='summer')
        #sns.heatmap(tabl['Percentage of Interaction'], annot=True)
        Highest_value= tabl['Percentage of Interaction'][tabl['Percentage of Interaction']==tabl['Percentage of Interaction'].max()]
        Highest_value=Highest_value.to_dict()
        print "Highest percenrage of Interactions identified","<br/>"
        Max_tabl=pd.DataFrame(Highest_value.items())
        Max_tabl.columns = ['Ligand Atom', 'Percentage']
        Max_tabl.rename(index={0: 'Highest'})
        #Max_tabl=pd.Series(Highest_value).to_frame()
        #Max_tabl.index.rename = 'index'
        #Max_tabl.rename(index={0:'zero'}, inplace=True)
        #df1.rename(index={0: 'a'})
        print Max_tabl.T.to_html(justify='center')
    else:
        print "No Interactions Observed"
######End of Percentage section###

####Start of Distance section##

#####Start of Distance section###
def distance_calc(dictnames):
    DistMean_dict={}
    DistFinal_pdb={}
    if bool(dictnames):
        for key,value in dictnames.iteritems():
            for key1,value1 in value.iteritems():
                results = map(float, value1)
                #print value1, np.mean(results)
                mean1=round(np.float64(np.mean(results)), 2)
                DistMean_dict[key1]=mean1
                DistFinal_pdb[key]=DistMean_dict
            DistMean_dict={}
        Distance_tabl=pd.DataFrame.from_dict(DistFinal_pdb)
        print Distance_tabl.T.to_html(justify='center'),"<br/>"
        print Distance_tabl.apply(pd.Series.describe,  axis=1)[['count','mean','std']].dropna().round(2).T.to_html(justify='center'),"<br/>"
        #Distance_tabl['Standard Deviation']=Distance_tabl.std(axis=1)
        #Distance_tabl['Standard Deviation']=Distance_tabl['Standard Deviation'].round(2)
        
        
    else:
        print "No Interactions Observed","<br/>"

#End of distance section##

####################End of Define function for Statistics ################################
aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

SubstructureExcel = str(uuid.uuid4())+'_CSV.csv'
writer = pd.ExcelWriter(SubstructureExcel, engine='xlsxwriter')

###Link to download file
print '<p style=text-align:center >Download: <a href=%s download>Interaction Data</a>'% SubstructureExcel

print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>Nicot sub group structure","</p>"
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. Compiled Bonded Interactions - Click to read Basic Statistics Information</button>"#Start of click drop down
print "<div class='contentsection'>"


print "<p style='font-size:20px; color:black' align='center'>"

print " Number of Ligand atoms:", len(Nicot), "<br/>"
print " Number of PDB IDs:", len(Nicot_allNH_Lig_Resdict.keys()), "<br/>"

print "<div class='row'>"# spliting into two columns
print "<div class='column'>"# spliting into two columns

if bool(Nicot_allH_Lig_Resdict):
    print "[1: Intercating and 0: Not Interacting]" , "<br/>"
    print percentage(Nicot_allH_Lig_Resdict, Nicot) 

if bool(Nicot_allH_Lig_Resdict_distance):
    print distance_calc(Nicot_allH_Lig_Resdict_distance) 

print "</div>"# closing of first columns

print "<div class='column'>"
if bool(Nicot_allNH_Lig_Resdict):
    print "Statistics of Non-Bonded Intercations"
    
    print  percentage(Nicot_allNH_Lig_Resdict,Nicot) 

if bool(Nicot_allNH_Lig_Resdict_distance):
    print distance_calc(Nicot_allNH_Lig_Resdict_distance) 

print "</div>"# closing of second columns
print "</div>"#closing of row

print "</div>"#End of click drop down
print "<br/>"

print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""

if bool(Nicot_allH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Nicot_allH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################All Residues Colored Table for Nicot: H bonded################################


####################All Residues Colored Table for Nicot: H bonded################################
H_templist4graph=[]
H_graphdic1={} 

if bool(Nicot_graphdicH):
    for k,v in Nicot_graphdicH.iteritems():
        #print k
        for value in v:
            H_templist4graph.append(value)
            samp=sorted(list(set(H_templist4graph)))
        H_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        H_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in H_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
    
    length_ofcell=max(length_listofcompiledresidues)   
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(H_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in H_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
if bool(Nicot_allNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print pd.DataFrame.from_dict(Nicot_allNH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print "No Interactions"
####################All Residues Colored Table for Nicot: NON bonded################################
NH_templist4graph=[]
NH_graphdic1={} 


if bool(Nicot_graphdicNH):
    for k,v in Nicot_graphdicNH.iteritems():
        #print k
        for value in v:
            NH_templist4graph.append(value)
            samp=sorted(list(set(NH_templist4graph)))
        NH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        #print temlist
        #print samp
        
        NH_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in NH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(NH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in NH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "No Interactions"


print """
        </div> 
  </div>  
"""#closing of col-2-3 and module


print """


   <div class="col-2-3">
     <div class="module">
   
"""# initializing the middle column
if bool(Nicot_CommonH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Nicot_CommonH_Lig_Resdict).to_html(justify='center')#for common ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Nicot : H bonded################################
CommH_templist4graph=[]
CommH_graphdic1={} 
if bool(Nicot_common_graphdicH):
    
    for k,v in Nicot_common_graphdicH.iteritems():
        for value in v:
            CommH_templist4graph.append(value)
            samp=sorted(list(set(CommH_templist4graph)))
        CommH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        CommH_templist4graph=[] 
    length_listofcompiled_Common_residues=[]
    for key,value in CommH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
        
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommH_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
    
if bool(Nicot_CommonNH_Lig_Resdict):
    
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print pd.DataFrame.from_dict(Nicot_CommonNH_Lig_Resdict).to_html(justify='center')#for Common ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Nicot: NON bonded################################
CommNH_templist4graph=[]
CommNH_graphdic1={} 
if bool(Nicot_common_graphdicNH):
    for k,v in Nicot_common_graphdicNH.iteritems():
        #print k
        for value in v:
            CommNH_templist4graph.append(value)
            samp=sorted(list(set(CommNH_templist4graph)))
        CommNH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        CommNH_templist4graph=[] 
    length_listofcompile_Common_dresidues=[]
    for key,value in CommNH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
    
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommNH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommNH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "No Interactions"
print """
        </div>
    </div>
"""# closinf of column and module div




###############Web logo for Common Residues Section: H bonding#######################
print """
   <div class="col-2-3">
     <div class="module">
     """
Nicot_graph_filename = str(uuid.uuid4())


Weblogo_dict_H={}
Weblogo_dict_H1={}
if bool (CommH_graphdic1):
    for key in sorted(CommH_graphdic1):
        for i in CommH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_H.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_H1.setdefault('%s'%m,{}).update(counted)
    
    
    
    zipfilename='/tmp/'+Nicot_graph_filename+'_Hbonding'+'.zip'
    
    Nicot_aminoacid_singlecode={}
    aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    for Nicot_ligand_key, Nicot_amino_frequency in Weblogo_dict_H1.iteritems():
        #print ligand_key
        for i in Nicot_ligand_key:
            
            for Nicot_amino,Nicot_frequency in Nicot_amino_frequency.iteritems():
                
                for Nicot_amino_3letter,Nicot_code_frequency in aminoacid_code.iteritems():
                    if Nicot_amino == Nicot_amino_3letter:
                        recoded[Nicot_code_frequency]=Nicot_frequency
                        Nicot_aminoacid_singlecode.setdefault('%s'%Nicot_ligand_key,{}).update(recoded)
                        recoded={}
    
    Nicot_Frequency=1  
    instances=[]
    Nicot_weblogo_collection=[]                  
    for Nicot_ligand_key1, amino_frequency1 in Nicot_aminoacid_singlecode.iteritems():
    
        for Nicot_Amino1, Nicot_number in amino_frequency1.iteritems():
    
            Nicot_Frequency=1
            while Nicot_Frequency <= Nicot_number:
                instances.append(Seq(Nicot_Amino1, IUPAC.protein))
                Nicot_Frequency=Nicot_Frequency+1
    
        Nicot_motif = motifs.create(instances)
    
        Nicot_mymotif ='/tmp/'+ Nicot_graph_filename+ '_H_'+ Nicot_ligand_key1 +'.svg'
        Nicot_motif.weblogo('%s'%Nicot_mymotif,format='SVG',xaxis_label= '%s' %Nicot_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Nicot_weblogo_collection.append(Nicot_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Nicot_weblogo_collection)
    


    
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    
    for Nicot_image in sorted(Nicot_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %Nicot_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Nicot_image
        print "</div>"
    
    print "</div>"
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Nicot_myzip:
        for Nicot_Images in Nicot_weblogo_collection:
            
            Nicot_myzip.write(Nicot_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Bonded Interactions:</p>"
    print "No Interactions"
###############Web logo for Common Residues Section: NON bonding#######################


Weblogo_dict_NH={}
Weblogo_dict_NH1={}
if bool(CommNH_graphdic1):
    for key in sorted(CommNH_graphdic1):
        for i in CommNH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_NH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_NH.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_NH1.setdefault('%s'%m,{}).update(counted)
    
    zipfilename='/tmp/'+Nicot_graph_filename+'_NHbonding'+'.zip'
    
    Nicot_aminoacid_singlecode={}
    
    recoded={}
    for Nicot_ligand_key, Nicot_amino_frequency in Weblogo_dict_NH1.iteritems():
        #print ligand_key
        for i in Nicot_ligand_key:
            
            for Nicot_amino,Nicot_frequency in Nicot_amino_frequency.iteritems():
                
                for Nicot_amino_3letter,Nicot_code_frequency in aminoacid_code.iteritems():
                    if Nicot_amino == Nicot_amino_3letter:
                        recoded[Nicot_code_frequency]=Nicot_frequency
                        Nicot_aminoacid_singlecode.setdefault('%s'%Nicot_ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Nicot_Frequency=1  
    instances=[]
    Nicot_weblogo_collection=[]                  
    for Nicot_ligand_key1, amino_frequency1 in Nicot_aminoacid_singlecode.iteritems():
    
        for Nicot_Amino1, Nicot_number in amino_frequency1.iteritems():
    
            Nicot_Frequency=1
            while Nicot_Frequency <= Nicot_number:
                instances.append(Seq(Nicot_Amino1, IUPAC.protein))
                Nicot_Frequency=Nicot_Frequency+1
    
        Nicot_motif = motifs.create(instances)
    
        Nicot_mymotif ='/tmp/'+ Nicot_graph_filename+ '_NH_'+ Nicot_ligand_key1 +'.svg'
        Nicot_motif.weblogo('%s'%Nicot_mymotif,format='SVG',xaxis_label= '%s' %Nicot_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Nicot_weblogo_collection.append(Nicot_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Nicot_weblogo_collection)
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>" #initiation of weblog_row
    
    
    for Nicot_image in sorted(Nicot_weblogo_collection):
        print "<div class='weblogo_column'>" #initiation of weblog_column
        print "<embed src='%s#page=1&view=FitH ' />" %Nicot_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='200' border='0'></iframe>"%Nicot_image
        print "</div>"#closing of weblog_column
    
    print "</div>"#closing of weblog_row
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Nicot_myzip:
        for Nicot_Images in Nicot_weblogo_collection:
            
            Nicot_myzip.write(Nicot_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Nonbonded Interactions:</p>"
    print "No Interactions"




print """
        </div>
    </div>


</div>
""" # closing of Nicot section



print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>Ribitol sub group structure","</p>"
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. Compiled Bonded Interactions - Click to read Basic Statistics Information</button>"#Start of click drop down
print "<div class='contentsection'>"


print "<p style='font-size:20px; color:black' align='center'>"

print " Number of Ligand atoms:", len(Ribitol), "<br/>"
print " Number of PDB IDs:", len(Ribitol_allNH_Lig_Resdict.keys()), "<br/>"

print "<div class='row'>"# spliting into two columns
print "<div class='column'>"# spliting into two columns

if bool(Ribitol_allH_Lig_Resdict):
    print "Statistics of Bonded Intercations" 
    print percentage(Ribitol_allH_Lig_Resdict,Ribitol) 

if bool(Ribitol_allH_Lig_Resdict_distance):
    print distance_calc(Ribitol_allH_Lig_Resdict_distance) 

print "</div>"# closing of first columns

print "<div class='column'>"
if bool(Ribitol_allNH_Lig_Resdict):
    print "Statistics of Non-Bonded Intercations"
    print  percentage(Ribitol_allNH_Lig_Resdict,Ribitol) 

if bool(Ribitol_allNH_Lig_Resdict_distance):
    print distance_calc(Ribitol_allNH_Lig_Resdict_distance) 

print "</div>"# closing of second columns
print "</div>"#closing of row

print "</div>"#End of click drop down
print "<br/>"



print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""

if bool(Ribitol_allH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Ribitol_allH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################All Residues Colored Table for Ribitol: H bonded################################


####################All Residues Colored Table for Ribitol: H bonded################################
H_templist4graph=[]
H_graphdic1={} 

if bool(Ribitol_graphdicH):
    for k,v in Ribitol_graphdicH.iteritems():
        #print k
        for value in v:
            H_templist4graph.append(value)
            samp=sorted(list(set(H_templist4graph)))
        H_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        H_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in H_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
    
    length_ofcell=max(length_listofcompiledresidues)   
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(H_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in H_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
if bool(Ribitol_allNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print pd.DataFrame.from_dict(Ribitol_allNH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print "No Interactions"
####################All Residues Colored Table for Ribitol: NON bonded################################
NH_templist4graph=[]
NH_graphdic1={} 


if bool(Ribitol_graphdicNH):
    for k,v in Ribitol_graphdicNH.iteritems():
        #print k
        for value in v:
            NH_templist4graph.append(value)
            samp=sorted(list(set(NH_templist4graph)))
        NH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        #print temlist
        #print samp
        
        NH_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in NH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(NH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in NH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "No Interactions"


print """
        </div> 
  </div>  
"""#closing of col-2-3 and module


print """


   <div class="col-2-3">
     <div class="module">
   
"""# initializing the middle column
if bool(Ribitol_CommonH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Ribitol_CommonH_Lig_Resdict).to_html(justify='center')#for common ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Ribitol : H bonded################################
CommH_templist4graph=[]
CommH_graphdic1={} 
if bool(Ribitol_common_graphdicH):
    
    for k,v in Ribitol_common_graphdicH.iteritems():
        for value in v:
            CommH_templist4graph.append(value)
            samp=sorted(list(set(CommH_templist4graph)))
        CommH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        CommH_templist4graph=[] 
    length_listofcompiled_Common_residues=[]
    for key,value in CommH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
        
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommH_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
    
if bool(Ribitol_CommonNH_Lig_Resdict):
    
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print pd.DataFrame.from_dict(Ribitol_CommonNH_Lig_Resdict).to_html(justify='center')#for Common ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Ribitol: NON bonded################################
CommNH_templist4graph=[]
CommNH_graphdic1={} 
if bool(Ribitol_common_graphdicNH):
    for k,v in Ribitol_common_graphdicNH.iteritems():
        #print k
        for value in v:
            CommNH_templist4graph.append(value)
            samp=sorted(list(set(CommNH_templist4graph)))
        CommNH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        CommNH_templist4graph=[] 
    length_listofcompile_Common_dresidues=[]
    for key,value in CommNH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
    
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommNH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommNH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "No Interactions"
print """
        </div>
    </div>
"""# closinf of column and module div




###############Web logo for Common Residues Section: H bonding#######################
print """
   <div class="col-2-3">
     <div class="module">
     """
Ribitol_graph_filename = str(uuid.uuid4())


Weblogo_dict_H={}
Weblogo_dict_H1={}
if bool (CommH_graphdic1):
    for key in sorted(CommH_graphdic1):
        for i in CommH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_H.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_H1.setdefault('%s'%m,{}).update(counted)
    
    
    
    zipfilename='/tmp/'+Ribitol_graph_filename+'_Hbonding'+'.zip'
    
    Ribitol_aminoacid_singlecode={}
    aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    for Ribitol_ligand_key, Ribitol_amino_frequency in Weblogo_dict_H1.iteritems():
        #print ligand_key
        for i in Ribitol_ligand_key:
            
            for Ribitol_amino,Ribitol_frequency in Ribitol_amino_frequency.iteritems():
                
                for Ribitol_amino_3letter,Ribitol_code_frequency in aminoacid_code.iteritems():
                    if Ribitol_amino == Ribitol_amino_3letter:
                        recoded[Ribitol_code_frequency]=Ribitol_frequency
                        Ribitol_aminoacid_singlecode.setdefault('%s'%Ribitol_ligand_key,{}).update(recoded)
                        recoded={}
    
    Ribitol_Frequency=1  
    instances=[]
    Ribitol_weblogo_collection=[]                  
    for Ribitol_ligand_key1, amino_frequency1 in Ribitol_aminoacid_singlecode.iteritems():
    
        for Ribitol_Amino1, Ribitol_number in amino_frequency1.iteritems():
    
            Ribitol_Frequency=1
            while Ribitol_Frequency <= Ribitol_number:
                instances.append(Seq(Ribitol_Amino1, IUPAC.protein))
                Ribitol_Frequency=Ribitol_Frequency+1
    
        Ribitol_motif = motifs.create(instances)
    
        Ribitol_mymotif ='/tmp/'+ Ribitol_graph_filename+ '_H_'+ Ribitol_ligand_key1 +'.svg'
        Ribitol_motif.weblogo('%s'%Ribitol_mymotif,format='SVG',xaxis_label= '%s' %Ribitol_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Ribitol_weblogo_collection.append(Ribitol_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Ribitol_weblogo_collection)
    


    
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    
    for Ribitol_image in sorted(Ribitol_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %Ribitol_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Ribitol_image
        print "</div>"
    
    print "</div>"
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Ribitol_myzip:
        for Ribitol_Images in Ribitol_weblogo_collection:
            
            Ribitol_myzip.write(Ribitol_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Bonded Interactions:</p>"
    print "No Interactions"
###############Web logo for Common Residues Section: NON bonding#######################


Weblogo_dict_NH={}
Weblogo_dict_NH1={}
if bool(CommNH_graphdic1):
    for key in sorted(CommNH_graphdic1):
        for i in CommNH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_NH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_NH.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_NH1.setdefault('%s'%m,{}).update(counted)
    
    zipfilename='/tmp/'+Ribitol_graph_filename+'_NHbonding'+'.zip'
    
    Ribitol_aminoacid_singlecode={}
    
    recoded={}
    for Ribitol_ligand_key, Ribitol_amino_frequency in Weblogo_dict_NH1.iteritems():
        #print ligand_key
        for i in Ribitol_ligand_key:
            
            for Ribitol_amino,Ribitol_frequency in Ribitol_amino_frequency.iteritems():
                
                for Ribitol_amino_3letter,Ribitol_code_frequency in aminoacid_code.iteritems():
                    if Ribitol_amino == Ribitol_amino_3letter:
                        recoded[Ribitol_code_frequency]=Ribitol_frequency
                        Ribitol_aminoacid_singlecode.setdefault('%s'%Ribitol_ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Ribitol_Frequency=1  
    instances=[]
    Ribitol_weblogo_collection=[]                  
    for Ribitol_ligand_key1, amino_frequency1 in Ribitol_aminoacid_singlecode.iteritems():
    
        for Ribitol_Amino1, Ribitol_number in amino_frequency1.iteritems():
    
            Ribitol_Frequency=1
            while Ribitol_Frequency <= Ribitol_number:
                instances.append(Seq(Ribitol_Amino1, IUPAC.protein))
                Ribitol_Frequency=Ribitol_Frequency+1
    
        Ribitol_motif = motifs.create(instances)
    
        Ribitol_mymotif ='/tmp/'+ Ribitol_graph_filename+ '_NH_'+ Ribitol_ligand_key1 +'.svg'
        Ribitol_motif.weblogo('%s'%Ribitol_mymotif,format='SVG',xaxis_label= '%s' %Ribitol_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Ribitol_weblogo_collection.append(Ribitol_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Ribitol_weblogo_collection)
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>" #initiation of weblog_row
    
    
    for Ribitol_image in sorted(Ribitol_weblogo_collection):
        print "<div class='weblogo_column'>" #initiation of weblog_column
        print "<embed src='%s#page=1&view=FitH ' />" %Ribitol_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='200' border='0'></iframe>"%Ribitol_image
        print "</div>"#closing of weblog_column
    
    print "</div>"#closing of weblog_row
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Ribitol_myzip:
        for Ribitol_Images in Ribitol_weblogo_collection:
            
            Ribitol_myzip.write(Ribitol_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Nonbonded Interactions:</p>"
    print "No Interactions"




print """
        </div>
    </div>


</div>
""" # closing of Ribitol section


print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>phosphate sub group structure","</p>"
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. Compiled Bonded Interactions - Click to read Basic Statistics Information</button>"#Start of click drop down
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"

print " Number of Ligand atoms:", len(phosphate), "<br/>"
print " Number of PDB IDs:", len(phosphate_allNH_Lig_Resdict.keys()), "<br/>"

print "<div class='row'>"# spliting into two columns
print "<div class='column'>"# spliting into two columns

if bool(phosphate_allH_Lig_Resdict):
    print "Statistics of Bonded Intercations" 
    print percentage(phosphate_allH_Lig_Resdict,phosphate) 

if bool(phosphate_allH_Lig_Resdict_distance):
    print distance_calc(phosphate_allH_Lig_Resdict_distance) 

print "</div>"# closing of first columns

print "<div class='column'>"

if bool(phosphate_allNH_Lig_Resdict):
    print "Statistics of Non-Bonded Intercations"
    print  percentage(phosphate_allNH_Lig_Resdict,phosphate) 

if bool(phosphate_allNH_Lig_Resdict_distance):
    print distance_calc(phosphate_allNH_Lig_Resdict_distance) 

print "</div>"# closing of second columns
print "</div>"#closing of row

print "</div>"#End of click drop down
print "<br/>"


print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""

if bool(phosphate_allH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(phosphate_allH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################All Residues Colored Table for phosphate: H bonded################################


####################All Residues Colored Table for phosphate: H bonded################################
H_templist4graph=[]
H_graphdic1={} 

if bool(phosphate_graphdicH):
    for k,v in phosphate_graphdicH.iteritems():
        #print k
        for value in v:
            H_templist4graph.append(value)
            samp=sorted(list(set(H_templist4graph)))
        H_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        H_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in H_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
    
    length_ofcell=max(length_listofcompiledresidues)   
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(H_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in H_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
if bool(phosphate_allNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print pd.DataFrame.from_dict(phosphate_allNH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print "No Interactions"
####################All Residues Colored Table for phosphate: NON bonded################################
NH_templist4graph=[]
NH_graphdic1={} 


if bool(phosphate_graphdicNH):
    for k,v in phosphate_graphdicNH.iteritems():
        #print k
        for value in v:
            NH_templist4graph.append(value)
            samp=sorted(list(set(NH_templist4graph)))
        NH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        #print temlist
        #print samp
        
        NH_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in NH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(NH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in NH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "No Interactions"


print """
        </div> 
  </div>  
"""#closing of col-2-3 and module


print """


   <div class="col-2-3">
     <div class="module">
   
"""# initializing the middle column
if bool(phosphate_CommonH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(phosphate_CommonH_Lig_Resdict).to_html(justify='center')#for common ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################Common Residues Colored Table for phosphate : H bonded################################
CommH_templist4graph=[]
CommH_graphdic1={} 
if bool(phosphate_common_graphdicH):
    
    for k,v in phosphate_common_graphdicH.iteritems():
        for value in v:
            CommH_templist4graph.append(value)
            samp=sorted(list(set(CommH_templist4graph)))
        CommH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        CommH_templist4graph=[] 
    length_listofcompiled_Common_residues=[]
    for key,value in CommH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
        
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommH_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"
    
if bool(phosphate_CommonNH_Lig_Resdict):
    
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print pd.DataFrame.from_dict(phosphate_CommonNH_Lig_Resdict).to_html(justify='center')#for Common ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print "No Interactions"
####################Common Residues Colored Table for phosphate: NON bonded################################
CommNH_templist4graph=[]
CommNH_graphdic1={} 
if bool(phosphate_common_graphdicNH):
    for k,v in phosphate_common_graphdicNH.iteritems():
        #print k
        for value in v:
            CommNH_templist4graph.append(value)
            samp=sorted(list(set(CommNH_templist4graph)))
        CommNH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        CommNH_templist4graph=[] 
    length_listofcompile_Common_dresidues=[]
    for key,value in CommNH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
    
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommNH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommNH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "No Interactions"
print """
        </div>
    </div>
"""# closinf of column and module div




###############Web logo for Common Residues Section: H bonding#######################
print """
   <div class="col-2-3">
     <div class="module">
     """
phosphate_graph_filename = str(uuid.uuid4())


Weblogo_dict_H={}
Weblogo_dict_H1={}
if bool (CommH_graphdic1):
    for key in sorted(CommH_graphdic1):
        for i in CommH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_H.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_H1.setdefault('%s'%m,{}).update(counted)
    
    
    
    zipfilename='/tmp/'+phosphate_graph_filename+'_Hbonding'+'.zip'
    
    phosphate_aminoacid_singlecode={}
    aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    for phosphate_ligand_key, phosphate_amino_frequency in Weblogo_dict_H1.iteritems():
        #print ligand_key
        for i in phosphate_ligand_key:
            
            for phosphate_amino,phosphate_frequency in phosphate_amino_frequency.iteritems():
                
                for phosphate_amino_3letter,phosphate_code_frequency in aminoacid_code.iteritems():
                    if phosphate_amino == phosphate_amino_3letter:
                        recoded[phosphate_code_frequency]=phosphate_frequency
                        phosphate_aminoacid_singlecode.setdefault('%s'%phosphate_ligand_key,{}).update(recoded)
                        recoded={}
    
    phosphate_Frequency=1  
    instances=[]
    phosphate_weblogo_collection=[]                  
    for phosphate_ligand_key1, amino_frequency1 in phosphate_aminoacid_singlecode.iteritems():
    
        for phosphate_Amino1, phosphate_number in amino_frequency1.iteritems():
    
            phosphate_Frequency=1
            while phosphate_Frequency <= phosphate_number:
                instances.append(Seq(phosphate_Amino1, IUPAC.protein))
                phosphate_Frequency=phosphate_Frequency+1
    
        phosphate_motif = motifs.create(instances)
    
        phosphate_mymotif ='/tmp/'+ phosphate_graph_filename+ '_H_'+ phosphate_ligand_key1 +'.svg'
        phosphate_motif.weblogo('%s'%phosphate_mymotif,format='SVG',xaxis_label= '%s' %phosphate_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        phosphate_weblogo_collection.append(phosphate_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in phosphate_weblogo_collection)
    


    
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    
    for phosphate_image in sorted(phosphate_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %phosphate_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%phosphate_image
        print "</div>"
    
    print "</div>"
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as phosphate_myzip:
        for phosphate_Images in phosphate_weblogo_collection:
            
            phosphate_myzip.write(phosphate_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Bonded Interactions:</p>"
    print "No Interactions"
###############Web logo for Common Residues Section: NON bonding#######################


Weblogo_dict_NH={}
Weblogo_dict_NH1={}
if bool(CommNH_graphdic1):
    for key in sorted(CommNH_graphdic1):
        for i in CommNH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_NH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_NH.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_NH1.setdefault('%s'%m,{}).update(counted)
    
    zipfilename='/tmp/'+phosphate_graph_filename+'_NHbonding'+'.zip'
    
    phosphate_aminoacid_singlecode={}
    
    recoded={}
    for phosphate_ligand_key, phosphate_amino_frequency in Weblogo_dict_NH1.iteritems():
        #print ligand_key
        for i in phosphate_ligand_key:
            
            for phosphate_amino,phosphate_frequency in phosphate_amino_frequency.iteritems():
                
                for phosphate_amino_3letter,phosphate_code_frequency in aminoacid_code.iteritems():
                    if phosphate_amino == phosphate_amino_3letter:
                        recoded[phosphate_code_frequency]=phosphate_frequency
                        phosphate_aminoacid_singlecode.setdefault('%s'%phosphate_ligand_key,{}).update(recoded)
                        recoded={}
    
    
    phosphate_Frequency=1  
    instances=[]
    phosphate_weblogo_collection=[]                  
    for phosphate_ligand_key1, amino_frequency1 in phosphate_aminoacid_singlecode.iteritems():
    
        for phosphate_Amino1, phosphate_number in amino_frequency1.iteritems():
    
            phosphate_Frequency=1
            while phosphate_Frequency <= phosphate_number:
                instances.append(Seq(phosphate_Amino1, IUPAC.protein))
                phosphate_Frequency=phosphate_Frequency+1
    
        phosphate_motif = motifs.create(instances)
    
        phosphate_mymotif ='/tmp/'+ phosphate_graph_filename+ '_NH_'+ phosphate_ligand_key1 +'.svg'
        phosphate_motif.weblogo('%s'%phosphate_mymotif,format='SVG',xaxis_label= '%s' %phosphate_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        phosphate_weblogo_collection.append(phosphate_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in phosphate_weblogo_collection)
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>" #initiation of weblog_row
    
    
    for phosphate_image in sorted(phosphate_weblogo_collection):
        print "<div class='weblogo_column'>" #initiation of weblog_column
        print "<embed src='%s#page=1&view=FitH ' />" %phosphate_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='200' border='0'></iframe>"%phosphate_image
        print "</div>"#closing of weblog_column
    
    print "</div>"#closing of weblog_row
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as phosphate_myzip:
        for phosphate_Images in phosphate_weblogo_collection:
            
            phosphate_myzip.write(phosphate_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo showing Common Nonbonded Interactions:</p>"
    print "No Interactions"




print """
        </div>
    </div>


</div>
""" # closing of phosphate section





#####################################################

print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>Ribose sub group structure","</p>"
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. Compiled Bonded Interactions - Click to read Basic Statistics Information</button>"#Start of click drop down
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"

print " Number of Ligand atoms:", len(Ribose), "<br/>"
print " Number of PDB IDs:", len(Ribose_allNH_Lig_Resdict.keys()), "<br/>"

print "<div class='row'>"# spliting into two columns
print "<div class='column'>"# spliting into two columns

if bool(Ribose_allH_Lig_Resdict):
    print "Statistics of Bonded Intercations" 
    print percentage(Ribose_allH_Lig_Resdict,Ribose) 

if bool(Ribose_allH_Lig_Resdict_distance):
    print distance_calc(Ribose_allH_Lig_Resdict_distance) 

print "</div>"# closing of first columns
print "<div class='column'>"

if bool(Ribose_allNH_Lig_Resdict):
    print "Statistics of Non-Bonded Intercations"
    
    print  percentage(Ribose_allNH_Lig_Resdict,Ribose) 

if bool(Ribose_allNH_Lig_Resdict_distance):
    print distance_calc(Ribose_allNH_Lig_Resdict_distance) 

print "</div>"# closing of second columns
print "</div>"#closing of row

print "</div>"#End of click drop down
print "<br/>"



print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
"""#start of    Ribose grid section

print Ribose_allH_Lig_Resdict
if bool(Ribose_allH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Ribose_allH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>"
    print "No Interactions"
####################All Residues Colored Table for Ribose: H bonded################################
H_templist4graph=[]
H_graphdic1={} 
if bool(Ribose_graphdicH):
    for k,v in Ribose_graphdicH.iteritems():
        #print k
        for value in v:
            H_templist4graph.append(value)
            samp=sorted(list(set(H_templist4graph)))
        H_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        H_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in H_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)   
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(H_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in H_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"


if bool(Ribose_allNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print pd.DataFrame.from_dict(Ribose_allNH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print "NO Interactions"
####################All Residues Colored Table for NON bonded################################
NH_templist4graph=[]
NH_graphdic1={} 
if bool(Ribose_graphdicNH):
    for k,v in Ribose_graphdicNH.iteritems():
        #print k
        for value in v:
            NH_templist4graph.append(value)
            samp=sorted(list(set(NH_templist4graph)))
        NH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        #print temlist
        #print samp
        
        NH_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in NH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(NH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in NH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "No Interactions"


print """
        </div> 
  </div>  
"""#closing of first col-2-3 and module


print """


   <div class="col-2-3">
     <div class="module">
   
""" #initializing of second column


if bool(Ribose_CommonH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Ribose_CommonH_Lig_Resdict).to_html(justify='center')#for common ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Ribose : H bonded################################


CommH_templist4graph=[]
CommH_graphdic1={}
if bool(Ribose_common_graphdicH): 
    for k,v in Ribose_common_graphdicH.iteritems():
        for value in v:
            CommH_templist4graph.append(value)
            samp=sorted(list(set(CommH_templist4graph)))
        CommH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        CommH_templist4graph=[] 
    length_listofcompiled_Common_residues=[]
    for key,value in CommH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
        
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommH_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"


if bool(Ribose_CommonNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print pd.DataFrame.from_dict(Ribose_CommonNH_Lig_Resdict).to_html(justify='center')#for Common ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Ribose: NON bonded################################


CommNH_templist4graph=[]
CommNH_graphdic1={} 
if bool(Ribose_common_graphdicNH):
    for k,v in Ribose_common_graphdicNH.iteritems():
        #print k
        for value in v:
            CommNH_templist4graph.append(value)
            samp=sorted(list(set(CommNH_templist4graph)))
        CommNH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        CommNH_templist4graph=[] 
    length_listofcompile_Common_dresidues=[]
    for key,value in CommNH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
    
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommNH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommNH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "No Interactions"
    
print """
        </div>
    </div>
"""# closinf of second column and module divi




###############Web logo for Common Residues Section: H bonding#######################
print """
   <div class="col-2-3">
     <div class="module">
     """
Ribose_graph_filename = str(uuid.uuid4())
Weblogo_dict_H={}
Weblogo_dict_H1={}
if bool (CommH_graphdic1):
    for key in sorted(CommH_graphdic1):
        for i in CommH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_H.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_H1.setdefault('%s'%m,{}).update(counted)
    
    
    
    zipfilename='/tmp/'+Ribose_graph_filename+'_Hbonding'+'.zip'
    
    Ribose_aminoacid_singlecode={}
    aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    for Ribose_ligand_key, Ribose_amino_frequency in Weblogo_dict_H1.iteritems():
        #print ligand_key
        for i in Ribose_ligand_key:
            
            for Ribose_amino,Ribose_frequency in Ribose_amino_frequency.iteritems():
                
                for Ribose_amino_3letter,Ribose_code_frequency in aminoacid_code.iteritems():
                    if Ribose_amino == Ribose_amino_3letter:
                        recoded[Ribose_code_frequency]=Ribose_frequency
                        Ribose_aminoacid_singlecode.setdefault('%s'%Ribose_ligand_key,{}).update(recoded)
                        recoded={}
    
    Ribose_Frequency=1  
    instances=[]
    Ribose_weblogo_collection=[]                  
    for Ribose_ligand_key1, amino_frequency1 in Ribose_aminoacid_singlecode.iteritems():
    
        for Ribose_Amino1, Ribose_number in amino_frequency1.iteritems():
    
            Ribose_Frequency=1
            while Ribose_Frequency <= Ribose_number:
                instances.append(Seq(Ribose_Amino1, IUPAC.protein))
                Ribose_Frequency=Ribose_Frequency+1
    
        Ribose_motif = motifs.create(instances)
    
        Ribose_mymotif ='/tmp/'+ Ribose_graph_filename+ '_H_'+ Ribose_ligand_key1 +'.svg'
        Ribose_motif.weblogo('%s'%Ribose_mymotif,format='SVG',xaxis_label= '%s' %Ribose_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Ribose_weblogo_collection.append(Ribose_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Ribose_weblogo_collection)
    


    
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    
    for Ribose_image in sorted(Ribose_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %Ribose_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Ribose_image
        print "</div>"
    
    print "</div>"
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Ribose_myzip:
        for Ribose_Images in Ribose_weblogo_collection:
            
            Ribose_myzip.write(Ribose_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo for Common Bonded Interactions:</p>"
    print "NO Interactions"
###############Web logo for Common Residues Section: NON bonding#######################


Weblogo_dict_NH={}
Weblogo_dict_NH1={}
if bool(CommNH_graphdic1):
    for key in sorted(CommNH_graphdic1):
        for i in CommNH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_NH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_NH.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_NH1.setdefault('%s'%m,{}).update(counted)
    
    zipfilename='/tmp/'+Ribose_graph_filename+'_NHbonding'+'.zip'
    
    Ribose_aminoacid_singlecode={}
    
    recoded={}
    for Ribose_ligand_key, Ribose_amino_frequency in Weblogo_dict_NH1.iteritems():
        #print ligand_key
        for i in Ribose_ligand_key:
            
            for Ribose_amino,Ribose_frequency in Ribose_amino_frequency.iteritems():
                
                for Ribose_amino_3letter,Ribose_code_frequency in aminoacid_code.iteritems():
                    if Ribose_amino == Ribose_amino_3letter:
                        recoded[Ribose_code_frequency]=Ribose_frequency
                        Ribose_aminoacid_singlecode.setdefault('%s'%Ribose_ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Ribose_Frequency=1  
    instances=[]
    Ribose_weblogo_collection=[]                  
    for Ribose_ligand_key1, amino_frequency1 in Ribose_aminoacid_singlecode.iteritems():
    
        for Ribose_Amino1, Ribose_number in amino_frequency1.iteritems():
    
            Ribose_Frequency=1
            while Ribose_Frequency <= Ribose_number:
                instances.append(Seq(Ribose_Amino1, IUPAC.protein))
                Ribose_Frequency=Ribose_Frequency+1
    
        Ribose_motif = motifs.create(instances)
    
        Ribose_mymotif ='/tmp/'+ Ribose_graph_filename+ '_NH_'+ Ribose_ligand_key1 +'.svg'
        Ribose_motif.weblogo('%s'%Ribose_mymotif,format='SVG',xaxis_label= '%s' %Ribose_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Ribose_weblogo_collection.append(Ribose_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Ribose_weblogo_collection)
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>" #initiation of weblog_row
    
    
    for Ribose_image in sorted(Ribose_weblogo_collection):
        print "<div class='weblogo_column'>" #initiation of weblog_column
        print "<embed src='%s#page=1&view=FitH ' />" %Ribose_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='200' border='0'></iframe>"%Ribose_image
        print "</div>"#closing of weblog_column
    
    print "</div>"#closing of weblog_row
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Ribose_myzip:
        for Ribose_Images in Ribose_weblogo_collection:
            
            Ribose_myzip.write(Ribose_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo for Common Nonbonded Interactions:</p>"
    print "No Interaction"




print """
        </div>
    </div>


</div>
""" # closing of RiboSE section

##############################

print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>Adenin sub group structure","</p>"
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. Compiled Bonded Interactions - Click to read Basic Statistics Information</button>"#Start of click drop down
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"

print " Number of Ligand atoms:", len(Adenin), "<br/>"
print " Number of PDB IDs:", len(Adenin_allNH_Lig_Resdict.keys()), "<br/>"

print "<div class='row'>"# spliting into two columns
print "<div class='column'>"# spliting into two columns

if bool(Adenin_allH_Lig_Resdict):
    print "Statistics of Bonded Intercations" 
    print percentage(Adenin_allH_Lig_Resdict,Adenin)
else:
    print "No data"

if bool(Adenin_allH_Lig_Resdict_distance):
    print distance_calc(Adenin_allH_Lig_Resdict_distance) 

print "</div>"# closing of first columns

print "<div class='column'>"

if bool(Adenin_allNH_Lig_Resdict):
    print "Statistics of Non-Bonded Intercations"
    
    print  percentage(Adenin_allNH_Lig_Resdict,Adenin)

if bool(Adenin_allNH_Lig_Resdict_distance):
    print distance_calc(Adenin_allNH_Lig_Resdict_distance) 

print "</div>"# closing of second columns
print "</div>"#closing of row

print "</div>"#End of click drop down
print "<br/>"


print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
"""#Initialization of    Adenin grid section


if bool(Adenin_allH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Adenin_allH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: hydrogen bonds contacts"  ,"</p>" 
    print "No Interactions"
####################All Residues Colored Table for Adenin: H bonded################################


H_templist4graph=[]
H_graphdic1={} 
if bool(Adenin_graphdicH):
    for k,v in Adenin_graphdicH.iteritems():
        #print k
        for value in v:
            H_templist4graph.append(value)
            samp=sorted(list(set(H_templist4graph)))
        H_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        H_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in H_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)   
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(H_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in H_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: hydrogen bonds contacts ","</p>"
    print "No Interactions"


if bool(Adenin_allNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print pd.DataFrame.from_dict(Adenin_allNH_Lig_Resdict).to_html(justify='center')#for all ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of residues: non-bonded contacts","</p>"
    print "No Interactions"
####################All Residues Colored Table for NON bonded################################
NH_templist4graph=[]
NH_graphdic1={} 
if bool(Adenin_graphdicNH):
    for k,v in Adenin_graphdicNH.iteritems():
        #print k
        for value in v:
            NH_templist4graph.append(value)
            samp=sorted(list(set(NH_templist4graph)))
        NH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        #print temlist
        #print samp
        
        NH_templist4graph=[] 
    length_listofcompiledresidues=[]
    for key,value in NH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            #print valu
            #print len(valu)
            length_listofcompiledresidues.append(len(valu))
            
    length_ofcell=max(length_listofcompiledresidues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(NH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in NH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of amino acids: non-bonded contacts","</p>"
    print "No Interactions"


print """
        </div> 
  </div>  
"""#closing of col-2-3 and module


print """


   <div class="col-2-3">
     <div class="module">
   
"""


if bool(Adenin_CommonH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print pd.DataFrame.from_dict(Adenin_CommonH_Lig_Resdict).to_html(justify='center')#for common ligand atoms - hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: hydrogen bonds contacts"  ,"</p>" 
    print "<p> No Common Interactions</p>"     
####################Common Residues Colored Table for Adenin : H bonded################################
CommH_templist4graph=[]
CommH_graphdic1={} 


if bool(Adenin_common_graphdicH):
    for k,v in Adenin_common_graphdicH.iteritems():
        for value in v:
            CommH_templist4graph.append(value)
            samp=sorted(list(set(CommH_templist4graph)))
        CommH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
        
        CommH_templist4graph=[] 
    length_listofcompiled_Common_residues=[]
    for key,value in CommH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
        
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    #print  "<br/>"
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommH_graphdic1[key]:
            dat1= g1.split(', ')
            for H_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if H_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%H_k3                    
                
                if H_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%H_k3
                    
                if H_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%H_k3
    
                if H_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%H_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"


else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based color-coding of common amino acids: hydrogen bonds contacts ","</p>"
    print "<p> No Common Atoms Identified</p>" 


if bool(Adenin_CommonNH_Lig_Resdict):
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print pd.DataFrame.from_dict(Adenin_CommonNH_Lig_Resdict).to_html(justify='center')#for Common ligand atoms - Non hydrogen bonded
else:
    print "<p style='font-size:20px; color:brown'>List of common residues: non-bonded contacts","</p>" 
    print "No Interactions"
####################Common Residues Colored Table for Adenin: NON bonded################################


CommNH_templist4graph=[]
CommNH_graphdic1={} 
if bool(Adenin_common_graphdicNH):
    for k,v in Adenin_common_graphdicNH.iteritems():
        #print k
        for value in v:
            CommNH_templist4graph.append(value)
            samp=sorted(list(set(CommNH_templist4graph)))
        CommNH_graphdic1.setdefault('%s'%k,[]).append(', '.join(samp))
    
        CommNH_templist4graph=[] 
    length_listofcompile_Common_dresidues=[]
    for key,value in CommNH_graphdic1.iteritems():    
        for i in value:
            valu=i.split(', ')
            length_listofcompiled_Common_residues.append(len(valu))
    
    
    length_ofcell=max(length_listofcompiled_Common_residues)
    
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>Compiled List of Common Residues From Analysed Protein Structures</th>"% length_ofcell
    print "</tr>"
    for key in sorted(CommNH_graphdic1.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1 in CommNH_graphdic1[key]:
            dat1= g1.split(', ')
            for NH_k3 in dat1:
                print "<td align='center'>"
                #print k3
                if NH_k3.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%NH_k3                    
                
                if NH_k3.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%NH_k3
                    
                if NH_k3.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%NH_k3
    
                if NH_k3.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%NH_k3
                
                print "</td>"
        #print "<tr>"    
        print "</tr>"
    print "</table>"
else:
    print "<p style='font-size:20px; color:brown'> Physicochemical property based coloring of Common amino acids:Nonbonded Contacts","</p>"
    print "No Interactions"


print """
        </div>
    </div>
"""# closinf of column and module divi




###############Web logo for Common Residues Section: H bonding#######################
print """
   <div class="col-2-3">
     <div class="module">
     """
Adenin_graph_filename = str(uuid.uuid4())
Weblogo_dict_H={}
Weblogo_dict_H1={}


if bool(CommH_graphdic1):
    for key in sorted(CommH_graphdic1):
        for i in CommH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_H.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_H1.setdefault('%s'%m,{}).update(counted)
    
    
    
    zipfilename='/tmp/'+Adenin_graph_filename+'_Hbonding'+'.zip'
    
    Adenin_aminoacid_singlecode={}
    aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    for Adenin_ligand_key, Adenin_amino_frequency in Weblogo_dict_H1.iteritems():
        #print ligand_key
        for i in Adenin_ligand_key:
            
            for Adenin_amino,Adenin_frequency in Adenin_amino_frequency.iteritems():
                
                for Adenin_amino_3letter,Adenin_code_frequency in aminoacid_code.iteritems():
                    if Adenin_amino == Adenin_amino_3letter:
                        recoded[Adenin_code_frequency]=Adenin_frequency
                        Adenin_aminoacid_singlecode.setdefault('%s'%Adenin_ligand_key,{}).update(recoded)
                        recoded={}
    
    Adenin_Frequency=1  
    instances=[]
    Adenin_weblogo_collection=[]                  
    for Adenin_ligand_key1, amino_frequency1 in Adenin_aminoacid_singlecode.iteritems():
    
        for Adenin_Amino1, Adenin_number in amino_frequency1.iteritems():
    
            Adenin_Frequency=1
            while Adenin_Frequency <= Adenin_number:
                instances.append(Seq(Adenin_Amino1, IUPAC.protein))
                Adenin_Frequency=Adenin_Frequency+1
    
        Adenin_motif = motifs.create(instances)
    
        Adenin_mymotif ='/tmp/'+ Adenin_graph_filename+ '_H_'+ Adenin_ligand_key1 +'.svg'
        Adenin_motif.weblogo('%s'%Adenin_mymotif,format='SVG',xaxis_label= '%s' %Adenin_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Adenin_weblogo_collection.append(Adenin_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Adenin_weblogo_collection)
    


    
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:</p>"
    
    print "<div class='weblogo_row'>"
    
    
    for Adenin_image in sorted(Adenin_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %Adenin_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Adenin_image
        print "</div>"
    
    print "</div>"
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Adenin_myzip:
        for Adenin_Images in Adenin_weblogo_collection:
            
            Adenin_myzip.write(Adenin_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo for Common bonded Interactions:</p>"
    print "No Interactions"
###############Web logo for Common Residues Section: NON bonding#######################


Weblogo_dict_NH={}
Weblogo_dict_NH1={}
if bool(CommNH_graphdic1):
    for key in sorted(CommNH_graphdic1):
        for i in CommNH_graphdic1[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Weblogo_dict_NH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Weblogo_dict_NH.iteritems():
        counted=dict(Counter(n))
        Weblogo_dict_NH1.setdefault('%s'%m,{}).update(counted)
    
    zipfilename='/tmp/'+Adenin_graph_filename+'_NHbonding'+'.zip'
    
    Adenin_aminoacid_singlecode={}
    
    recoded={}
    for Adenin_ligand_key, Adenin_amino_frequency in Weblogo_dict_NH1.iteritems():
        #print ligand_key
        for i in Adenin_ligand_key:
            
            for Adenin_amino,Adenin_frequency in Adenin_amino_frequency.iteritems():
                
                for Adenin_amino_3letter,Adenin_code_frequency in aminoacid_code.iteritems():
                    if Adenin_amino == Adenin_amino_3letter:
                        recoded[Adenin_code_frequency]=Adenin_frequency
                        Adenin_aminoacid_singlecode.setdefault('%s'%Adenin_ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Adenin_Frequency=1  
    instances=[]
    Adenin_weblogo_collection=[]                  
    for Adenin_ligand_key1, amino_frequency1 in Adenin_aminoacid_singlecode.iteritems():
    
        for Adenin_Amino1, Adenin_number in amino_frequency1.iteritems():
    
            Adenin_Frequency=1
            while Adenin_Frequency <= Adenin_number:
                instances.append(Seq(Adenin_Amino1, IUPAC.protein))
                Adenin_Frequency=Adenin_Frequency+1
    
        Adenin_motif = motifs.create(instances)
    
        Adenin_mymotif ='/tmp/'+ Adenin_graph_filename+ '_NH_'+ Adenin_ligand_key1 +'.svg'
        Adenin_motif.weblogo('%s'%Adenin_mymotif,format='SVG',xaxis_label= '%s' %Adenin_ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Adenin_weblogo_collection.append(Adenin_mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Adenin_weblogo_collection)
    print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:</p>"
    
    print "<div class='weblogo_row'>" #initiation of weblog_row
    
    
    for Adenin_image in sorted(Adenin_weblogo_collection):
        print "<div class='weblogo_column'>" #initiation of weblog_column
        print "<embed src='%s#page=1&view=FitH ' />" %Adenin_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='200' border='0'></iframe>"%Adenin_image
        print "</div>"#closing of weblog_column
    
    print "</div>"#closing of weblog_row
    
    ####zip file
        
    with ZipFile('%s'%zipfilename, 'w') as Adenin_myzip:
        for Adenin_Images in Adenin_weblogo_collection:
            
            Adenin_myzip.write(Adenin_Images)
else:
    print "<p style='font-size:20px; color:brown'> Weblogo for Common Nonbonded Interactions:</p>"
    print "No Interactions"


###########To write the dataframes to excel for download
Nicot_allH=pd.DataFrame.from_dict(Nicot_allH_Lig_Resdict)
Nicot_allH.to_excel(writer, sheet_name='Nicot_allH')
Nicot_allNH=pd.DataFrame.from_dict(Nicot_allNH_Lig_Resdict)
Nicot_allNH.to_excel(writer, sheet_name='Nicot_allNH')
Nicot_CommonH=pd.DataFrame.from_dict(Nicot_CommonH_Lig_Resdict)
Nicot_CommonH.to_excel(writer, sheet_name='Nicot_CommonH')
Nicot_CommonNH=pd.DataFrame.from_dict(Nicot_CommonNH_Lig_Resdict)
Nicot_CommonNH.to_excel(writer, sheet_name='Nicot_CommonNH')
Ribitol_allH=pd.DataFrame.from_dict(Ribitol_allH_Lig_Resdict)
Ribitol_allH.to_excel(writer, sheet_name='Ribitol_allH')
Ribitol_allNH=pd.DataFrame.from_dict(Ribitol_allNH_Lig_Resdict)
Ribitol_allNH.to_excel(writer, sheet_name='Ribitol_allNH')
Ribitol_CommonH=pd.DataFrame.from_dict(Ribitol_CommonH_Lig_Resdict)
Ribitol_CommonH.to_excel(writer, sheet_name='Ribitol_CommonH')
Ribitol_CommonNH=pd.DataFrame.from_dict(Ribitol_CommonNH_Lig_Resdict)
Ribitol_CommonNH.to_excel(writer, sheet_name='Ribitol_CommonNH')
phosphate_allH=pd.DataFrame.from_dict(phosphate_allH_Lig_Resdict)
phosphate_allH.to_excel(writer, sheet_name='phosphate_allH')
phosphate_allNH=pd.DataFrame.from_dict(phosphate_allNH_Lig_Resdict)
phosphate_allNH.to_excel(writer, sheet_name='phosphate_allNH')
phosphate_CommonH=pd.DataFrame.from_dict(phosphate_CommonH_Lig_Resdict)
phosphate_CommonH.to_excel(writer, sheet_name='phosphate_CommonH')
phosphate_CommonNH=pd.DataFrame.from_dict(phosphate_CommonNH_Lig_Resdict)
phosphate_CommonNH.to_excel(writer, sheet_name='phosphate_CommonNH')
Ribose_allH=pd.DataFrame.from_dict(Ribose_allH_Lig_Resdict)
Ribose_allH.to_excel(writer, sheet_name='Ribose_allH')
Ribose_allNH=pd.DataFrame.from_dict(Ribose_allNH_Lig_Resdict)
Ribose_allNH.to_excel(writer, sheet_name='Ribose_allNH')
Ribose_CommonH=pd.DataFrame.from_dict(Ribose_CommonH_Lig_Resdict)
Ribose_CommonH.to_excel(writer, sheet_name='Ribose_CommonH')
Ribose_CommonNH=pd.DataFrame.from_dict(Ribose_CommonNH_Lig_Resdict)
Ribose_CommonNH.to_excel(writer, sheet_name='Ribose_CommonNH')
Adenin_allH=pd.DataFrame.from_dict(Adenin_allH_Lig_Resdict)
Adenin_allH.to_excel(writer, sheet_name='Adenin_allH')
Adenin_allNH=pd.DataFrame.from_dict(Adenin_allNH_Lig_Resdict)
Adenin_allNH.to_excel(writer, sheet_name='Adenin_allNH')
Adenin_CommonH=pd.DataFrame.from_dict(Adenin_CommonH_Lig_Resdict)
Adenin_CommonH.to_excel(writer, sheet_name='Adenin_CommonH')
Adenin_CommonNH=pd.DataFrame.from_dict(Adenin_CommonNH_Lig_Resdict)
Adenin_CommonNH.to_excel(writer, sheet_name='Adenin_CommonNH')

writer.save()

print """
        </div>
    </div>


</div>
""" # closing of Adenin section

####Java script####
print """
<script>
var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
}
</script>
"""
###################   

print "</body>"
print "</html>"    
