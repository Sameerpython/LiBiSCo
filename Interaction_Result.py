#!/usr/bin/python

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
from collections import OrderedDict
#import seaborn as sns


from Bio.Alphabet import IUPAC
from zipfile import ZipFile
matplotlib.use('Agg')
pd.set_option('display.max_colwidth', -1)

# Create instance of FieldStorage
form = cgi.FieldStorage()


print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"
print "<style>"
# menubar style
print "ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}"
print "li{float:left;}"
print "li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}"
print "li a:hover { background-color: #111111;}"

# Grids style
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
    width: 20%;
    padding: 2px;
     overflow: scroll;
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
print "</style>"

print "<title>CoFactComp</title>"

###JMOLscript source
print "<script type='text/javascript' src='Jmol.js'></script>"
print "<script type='text/javascript' src='http://d3js.org/d3.v3.min.js'></script>"
print "<script type='text/javascript' src='http://mpld3.github.io/js/mpld3.v0.2.js'></script>"
#####End of JMOLscript source
print "</head>"
print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
print "<div align='center'>"
print "<img src='Title_image1.png' align='middle' width='1000' height='300'"
print "</div>"
print "<body>"
print "<ul>"
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href='Contact.py'>Contact</a></li>"
print "</ul>"
print "<div align='center'>"
print "<h2> Interaction Information </h2>" 
print "</div>"


#if form.getvalue('ligatom'):
#    lig_content = form.getvalue('ligatom')
#Ligscompare=list(lig_content.split(','))

# Information of the selected ligands and PDB ids from LigPage.py
variable = ""
value = ""
r = ""
value_dict={}
lig_sel=[]
for key in form.keys():
    if key !='ligatom':
        variable = str(key)
        value = str(form.getvalue(variable))
        value_dict.setdefault('%s'%variable,[]).append(value)
        r += "<p>"+ variable +", "+ value +"</p>\n"
print "<p style='font-size:20px; color:blue'> Results for the selected PDBID's and Ligands: ",', '.join("{}:{}".format(k,v) for k,v in value_dict.items()),"</p>","<br/>"

#PDBSUM URLs for connection
pdbsum_URL="http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode="
pdbsum_URL2="&template=links.html"

#DIctionary and List
pdbsum_dict=OrderedDict()
PDBID_LIST=[]

#Title for Page

Title="The Results are for the following selected PDBID's and Ligands:"

#Preparing PdbSum Url with selected PDB ids
for id,lig in value_dict.iteritems():
        pdbsumurl=pdbsum_URL+id+pdbsum_URL2
        lig_sel.append(lig)
        pdbsum_dict.setdefault('%s'%id,[]).append(pdbsumurl)

#creating a list for PDB ids
for id,url in pdbsum_dict.iteritems():
        PDBID_LIST.append(id)
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
URLlink_append=[]
for id,url in pdbsum_dict.iteritems():
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

#Looping over Ligand and PDBSum Urls (from above step) to extract the PDBSUm URL for the seleted Ligand Page in PDBSUm. The Ligand Page URL is now as a SET data type
for lig in lig_sel:

    for y in new:
        lig= ''.join(lig)
        if y.startswith(lig):
            y=y.split()
            link=y[1]
            link1="http://"+link

            link_set.add(link1)
            URLlink_append.append(link1)


    link_setlist=list(link_set)
#Using Beautifulsoup to extract all the links from PDBSUM Ligand interaction Page for each of the PDB ids.

PDBID_URL_dict=zip(PDBID_LIST,URLlink_append)
LiginteractPage=dict(PDBID_URL_dict)
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

                        final=i.split()
                        final=''.join(final)
                        
                        finalLIG_link.append(final)
                        lastitem=finalLIG_link[-1]
                        lastitem="http://"+lastitem

                if lastitem not in ligintelist:
                        ligintelist.append(lastitem)
                liginte_set.add(lastitem)
        liginte_list=list(liginte_set)
PDBID_INTURL_dict=zip(PDBID_LIST,ligintelist)
pdbsum_dict1=dict(PDBID_INTURL_dict)
for id,link in pdbsum_dict1.iteritems():
        links=list(pdbsum_dict1.viewvalues())
        PDBID=list(pdbsum_dict1.viewkeys())
Number_of_Ids=len(PDBID)
#FInal DIctionary with PDBID and Ligplot URL for extracting intercation details
mydictcheck={}
for ids,links in zip(PDBID,ligintelist):
        mydictcheck.setdefault('%s'%ids,[]).append(links)

#End of FInal DIctionary with PDBID and Ligplot URL for extracting intercation details

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
#End of selecting common ligand atoms that are non-hydrogen bonded in selected PDB structures




#extracting common atom names for the selected PDB ids (This combines both htdrogen and non bonding features.
atoms_commoncomp={}
for pdbids,pdbidlinks in mydictcheck.iteritems():
    for links_sel in pdbidlinks:
        links_sel1=str(links_sel)
        weblink=requests.get(links_sel1, stream=True)
        for atomlines in weblink.iter_lines():
            atomlines1=atomlines.strip()
            if atomlines1.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                atomlines2=atomlines1.split()
                atm_sel=atomlines2[8]
                atoms_commoncomp.setdefault('%s'%pdbids,[]).append(atm_sel)

atomsvalues_dict1=atoms_commoncomp.values()
common_intersectionfinal=sorted(list(set.intersection(*map(set,atomsvalues_dict1))))

#picking up the PDBid with the longest list of atomes
Maxlengthy_key=max(atoms_commoncomp, key=lambda m:len(set(atoms_commoncomp[m])))
valueoflengthyKEY=atoms_commoncomp.get(Maxlengthy_key)
valuesoflengthyKEY=sorted(list(set(valueoflengthyKEY)))
#End of this section picking up the PDBid with the longest list of atomes


#Merging of the atoms list for each PDB id into a a single list
listofAllatoms=[]
for k2,v2 in atoms_commoncomp.iteritems():
    for l2 in v2:
        atomitems=''.join(l2)
        listofAllatoms.append(atomitems)
finallistofAllatoms=sorted(list(set(listofAllatoms)))
#End of Merging of the atoms list for each PDB id into a a single list


listprint=[]
listprintdict={}
appended_lig_tabledic={}
appended_lig_tabledic1={}
if len(common_intersectionfinal)==0:
    print "<h1 style='color:blue'>No common atoms of the selected ligands are found to be interacting for the selected PDB ids.</h1>" 
    print "<h1 style='color:blue'>You can rerun the analysis by unselecting certain PDB ids.  </h1>"
    print "<h1 style='color:blue'>Below is the table listing the atoms of ligand each PDB id interacts.  </h1>"
    for k,v in atoms_commoncomp.iteritems():
        for i in k:
            pdb= k
            key= v
            for c1 in finallistofAllatoms:
                for b1 in v:
                    if c1 == b1:
                        if c1 not in listprint:
                            listprint.append(c1)
                            listprintdict.setdefault('%s'%pdb,[]).append(c1)
        listprint=[]
    df = pd.DataFrame.from_dict(listprintdict, orient='index').T.dropna()

    print df.apply(lambda c: pd.Series(c.name, c.values)).fillna('-').T.to_html(justify='center')

    print "<br/>"
    
    ##################Option for User to reselect PDB ids based on above table
    print "<h1 style='color:blue'>Please reselect PDB ids based on above observation for redoing the analysis.  </h1>"
    print "<form action='Interaction_reselection.py' method = 'post' target = '_blank'>"
    print "<table style=width:50%>"
    print "<tr>"
    print "<th>PDB ID</th>"
    print "<th colspan=2>LIGANDS</th>"
    for k,dk in value_dict.iteritems():
        count= len(dk)
        print "<tr>"
        print "<th rowspan='%d'>"% count,k,"</th>"
        for x in dk:
            print "<td align=center>", "%s" % x,"</td>"
            print "<td align=center>"
    	    print  "<input type='checkbox' name='%s' value='%s'/>" % (k,x)
    	    print "</td>"
    	    print "</tr>"

    print "</table>"
    print "<br/>"
    print "<input type = 'submit'  value = 'Submit'  />"
    print "<input type = 'reset'  value = 'Clear'  />"
    print "</form>"
            
else:
#hbond lists, dictionaries, sets initialised
    H_lresidue=[]
    H_latom=[]
    H_ldistance=[]
    residue={}
    atmname={}
    
    H_finalset=set()
    
    H_appended_lig_tabledic={}
    H_appended_lig_tabledic_withoutNAN={}
    printing = False
    H_combines_listdata=[]
    H_combines_listdata_withoutNAN=[]
    common_graphdic_H={}
    All_graphdic_H={}
    
    All_combine_Lig_Res_H={}
    AllH_combine_Lig_Res_distance={}
    AllH_combine_Lig_Res_distance_uniquify={}
    All_combine_Lig_Res_H_uniquify={}
    AllH_Lig_Resdict={}
    AllH_distance_Lig_dict={}

    Common_combine_Lig_Res_H={}
    Common_combine_Lig_Res_H_uniquify={}
    CommonH_Lig_Resdict={}
    
    CommonH_combine_Lig_Res_distance={}
    CommonH_combine_Lig_Res_distance_uniquify={}
    CommonH_distance_Lig_dict={}
    #selected ligandatoms
    LigAtom_optiongraphdic_H={}
    LigAtom_option_Lig_Res_H={}
    LigAtom_option_H_uniquify={}
    LigAtom_option_H_Lig_Resdict={}
    
    LigAtom_option_H_Lig_Res_distance={}
    LigAtom_option_H_Lig_Res_distance_uniquify={}
    LigAtom_option_H_distance_Lig_dict={}
    
    LigAtomOp_common_graphdic_H={}
    LigAtomOp_Common_combine_Lig_Res_H={}
    LigAtomOp_Common_combine_Lig_Res_H_uniquify={}
    LigAtomOp_CommonH_Lig_Resdict={}
    
    LigAtomOp_CommonH_combine_Lig_Res_distance={}
    LigAtomOp_CommonH_combine_Lig_Res_distance_uniquify={}
    LigAtomOp_CommonH_distance_Lig_dict={}


#Nonbonded lists, dictionaries, sets initialised
    NonH_lresidue=[]
    NonH_latom=[]
    NonH_ldistance=[]
    residue={}
    atmname={}
    
    NonH_finalset=set()
    
    NonH_appended_lig_tabledic={}
    NonH_appended_lig_tabledic_withoutNAN={}
    Nonh_printing = False
    NonH_combines_listdata_N=[]
    NonH_combines_listdata_withoutNAN=[] #change
    common_graphdic_NonH={}
    All_graphdic_NonH={}

    All_combine_Lig_Res_NH={}
    All_combine_Lig_Res_NH_uniquify={}
    AllNH_Lig_Resdict={}
    AllNH_combine_Lig_Res_distance={}
    AllNH_combine_Lig_Res_distance_uniquify={}
    AllNH_distance_Lig_dict={}



    Common_combine_Lig_Res_NH={}
    Common_combine_Lig_Res_NH_uniquify={}
    CommonNH_Lig_Resdict={}

    CommonNH_combine_Lig_Res_distance={}
    CommonNH_combine_Lig_Res_distance_uniquify={}
    CommonNH_distance_Lig_dict={}
    #selected ligandatoms
    LigAtom_optiongraphdic_NH={}
    LigAtom_option_Lig_Res_NH={}
    LigAtom_option_NH_uniquify={}
    LigAtom_option_NH_Lig_Resdict={}
    
    LigAtom_option_NH_Lig_Res_distance={}
    LigAtom_option_NH_Lig_Res_distance_uniquify={}
    LigAtom_option_NH_distance_Lig_dict={}
    
    LigAtomOp_common_graphdic_NH={}
    LigAtomOp_Common_combine_Lig_Res_NH={}
    LigAtomOp_Common_combine_Lig_Res_NH_uniquify={}
    LigAtomOp_CommonNH_Lig_Resdict={}
    
    LigAtomOp_CommonNH_combine_Lig_Res_distance={}
    LigAtomOp_CommonNH_combine_Lig_Res_distance_uniquify={}
    LigAtomOp_CommonNH_distance_Lig_dict={}


    for id,link in mydictcheck.iteritems():
    
       links_sel=link[0]
       link1= ''.join(str(links_sel))
       res2=urllib.urlopen(str(link1))
       html=res2.read()
    
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
                                line=line.split()
                                ligname=line[9]
                                atm=line[8]
                                res=line[3]
                                residuenum=line[4]
                                proteinchain=line[5]
                                distance=line[12]
                                resnum=res+residuenum+proteinchain
                                #appending each residue and its position to list called lresidue
                                H_lresidue.append(resnum)
                                #appending each ligand atom to list called latom
                                H_latom.append(atm)
                                #appending distance of each interaction to ldistance
                                H_ldistance.append(distance)
                                All_graphdic_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with all lig atom and residues for physio and weblogo
                                All_combine_Lig_Res_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with all lig atom and residues for table
                                All_combine_Lig_Res_H_uniquify= {k:list(set(j)) for k,j in All_combine_Lig_Res_H.items()}
                                AllH_Lig_Resdict['%s'%id]=All_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for all group
                                
                                AllH_combine_Lig_Res_distance.setdefault('%s'%atm,[]).append(distance)#creating dictionary with all lig atom and distance for table
                                AllH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in AllH_combine_Lig_Res_distance.items()}
                                AllH_distance_Lig_dict['%s'%id]=AllH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                if atm in H_common_intersectionfinal: # listing of common atms
                                    common_graphdic_H.setdefault('%s'%atm,[]).append(resnum)

                                    Common_combine_Lig_Res_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with common lig atom and residues for table
                                    Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in Common_combine_Lig_Res_H.items()}
                                    CommonH_Lig_Resdict['%s'%id]=Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group
                                    
                                    CommonH_combine_Lig_Res_distance.setdefault('%s'%atm,[]).append(distance)#creating dictionary with all lig atom and distance for table
                                    CommonH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in CommonH_combine_Lig_Res_distance.items()}
                                    CommonH_distance_Lig_dict['%s'%id]=CommonH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

 #                               if atm in Ligscompare:
 #                                   LigAtom_optiongraphdic_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with all lig atom and residues for physio and weblogo
 #                                   LigAtom_option_Lig_Res_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with all lig atom and residues for table
 #                                   LigAtom_option_H_uniquify= {k:list(set(j)) for k,j in LigAtom_option_Lig_Res_H.items()}
 #                                   LigAtom_option_H_Lig_Resdict['%s'%id]=LigAtom_option_H_uniquify#final dic for table with pdb id , lig atom and residues for all group

 #                                   LigAtom_option_H_Lig_Res_distance.setdefault('%s'%atm,[]).append(distance)#creating dictionary with all lig atom and distance for table
 #                                   LigAtom_option_H_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in LigAtom_option_H_Lig_Res_distance.items()}
 #                                   LigAtom_option_H_distance_Lig_dict['%s'%id]=LigAtom_option_H_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                    
 #                                   if atm in H_common_intersectionfinal:
 #                                      LigAtomOp_common_graphdic_H.setdefault('%s'%atm,[]).append(resnum)
 #                                      LigAtomOp_Common_combine_Lig_Res_H.setdefault('%s'%atm,[]).append(resnum)#creating dictionary with common lig atom and residues for table
 #                                      LigAtomOp_Common_combine_Lig_Res_H_uniquify={k:list(set(j)) for k,j in LigAtomOp_Common_combine_Lig_Res_H.items()}
 #                                      LigAtomOp_CommonH_Lig_Resdict['%s'%id]=LigAtomOp_Common_combine_Lig_Res_H_uniquify#final dic for table with pdb id , lig atom and residues for common group
                                       
 #                                      LigAtomOp_CommonH_combine_Lig_Res_distance.setdefault('%s'%atm,[]).append(distance)#creating dictionary with all lig atom and distance for table
 #                                      LigAtomOp_CommonH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in LigAtomOp_CommonH_combine_Lig_Res_distance.items()}
 #                                      LigAtomOp_CommonH_distance_Lig_dict['%s'%id]=LigAtomOp_CommonH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                       



                        else:
                            if line.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                                line=line.split()
                                ligname=line[9]
                                atmNH=line[8]
                                resNH=line[3]
                                residuenumNH=line[4]
                                proteinchainNH=line[5]
                                distanceNH=line[12]
                                resnumNH=resNH+residuenumNH+proteinchainNH
                                #appending each residue and its position to list called lresidue
                                NonH_lresidue.append(resnumNH)
                                #appending each ligand atom to list called latom
                                NonH_latom.append(atmNH)
                                #appending distance of each interaction to ldistance
                                NonH_ldistance.append(distance)
                                All_graphdic_NonH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for physio and weblogo

                                All_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
                                All_combine_Lig_Res_NH_uniquify= {k:list(set(j)) for k,j in All_combine_Lig_Res_NH.items()}
                                AllNH_Lig_Resdict['%s'%id]=All_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

                                AllNH_combine_Lig_Res_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and distance for table
                                AllNH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in AllNH_combine_Lig_Res_distance.items()}
                                AllNH_distance_Lig_dict['%s'%id]=AllNH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group


                                if atmNH in NONHcommon_intersectionfinal: # listing of common atms
                                    common_graphdic_NonH.setdefault('%s'%atmNH,[]).append(resnumNH)

                                    Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
                                    Common_combine_Lig_Res_NH_uniquify={k:list(set(j)) for k,j in Common_combine_Lig_Res_NH.items()}
                                    CommonNH_Lig_Resdict['%s'%id]=Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group

                                    CommonNH_combine_Lig_Res_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and distance for table
                                    CommonNH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in CommonNH_combine_Lig_Res_distance.items()}
                                    CommonNH_distance_Lig_dict['%s'%id]=CommonNH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group


#                                if atmNH in Ligscompare:
 #                                   LigAtom_optiongraphdic_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for physio and weblogo
 #                                   LigAtom_option_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with all lig atom and residues for table
 #                                   LigAtom_option_NH_uniquify= {k:list(set(j)) for k,j in LigAtom_option_Lig_Res_NH.items()}
 #                                   LigAtom_option_NH_Lig_Resdict['%s'%id]=LigAtom_option_NH_uniquify#final dic for table with pdb id , lig atom and residues for all group

 #                                   LigAtom_option_NH_Lig_Res_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and distance for table
 #                                   LigAtom_option_NH_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in LigAtom_option_NH_Lig_Res_distance.items()}
 #                                   LigAtom_option_NH_distance_Lig_dict['%s'%id]=LigAtom_option_NH_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group
                                    
 #                                   if atmNH in NONHcommon_intersectionfinal:
 #                                      LigAtomOp_common_graphdic_NH.setdefault('%s'%atmNH,[]).append(resnumNH)
 #                                      LigAtomOp_Common_combine_Lig_Res_NH.setdefault('%s'%atmNH,[]).append(resnumNH)#creating dictionary with common lig atom and residues for table
 #                                      LigAtomOp_Common_combine_Lig_Res_NH_uniquify={k:list(set(j)) for k,j in LigAtomOp_Common_combine_Lig_Res_NH.items()}
 #                                      LigAtomOp_CommonNH_Lig_Resdict['%s'%id]=LigAtomOp_Common_combine_Lig_Res_NH_uniquify#final dic for table with pdb id , lig atom and residues for common group
                                       
 #                                      LigAtomOp_CommonNH_combine_Lig_Res_distance.setdefault('%s'%atmNH,[]).append(distanceNH)#creating dictionary with all lig atom and distance for table
 #                                      LigAtomOp_CommonNH_combine_Lig_Res_distance_uniquify= {k:list(set(j)) for k,j in LigAtomOp_CommonNH_combine_Lig_Res_distance.items()}
 #                                      LigAtomOp_CommonNH_distance_Lig_dict['%s'%id]=LigAtomOp_CommonNH_combine_Lig_Res_distance_uniquify#final dic for table with pdb id , lig atom and distance for all group

#emptying hbond lists and sets

       H_lresidue=[]
       H_latom=[]
       H_ldistance=[]
       H_combines_listdata=[]
       H_combines_listdata_withoutNAN=[]
       H_finalset.clear()
       All_combine_Lig_Res_H={}
       Common_combine_Lig_Res_H={}
       AllH_combine_Lig_Res_distance={}
       CommonH_combine_Lig_Res_distance={}
       #ligand atom selected
       LigAtom_option_Lig_Res_H={}
       LigAtom_option_H_Lig_Res_distance={}
       LigAtomOp_Common_combine_Lig_Res_H={}
       LigAtomOp_CommonH_combine_Lig_Res_distance

#emptying Nonbonded lists and sets
       NonH_lresidue=[]
       NonH_latom=[]
       NonH_ldistance=[]
       NonH_combines_listdata_N=[]
       NonH_combines_listdata_withoutNAN=[]
       NonH_finalset.clear()
       All_combine_Lig_Res_NH={}
       Common_combine_Lig_Res_NH={}
       AllNH_combine_Lig_Res_distance={}
       CommonNH_combine_Lig_Res_distance={}
       #ligand atom selected
       LigAtom_option_Lig_Res_NH={}
       LigAtom_option_NH_Lig_Res_distance={}
       LigAtomOp_Common_combine_Lig_Res_NH={}
       LigAtomOp_CommonNH_combine_Lig_Res_distance

####################Define function for Statistics ################################
def percentage(dictname):
    
    Count_Atom={}
    percentage_Atom={}
    atmlist=[]
    if bool(dictname):
        for key, value in dictname.iteritems():
            
            for atom in finallistofAllatoms:
            
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
        print "<br/>"," No. of ligand atoms:", len(count_atmlist), "/",len(finallistofAllatoms), "<br/>"
        print tabl.T.to_html(justify='center'),"<br/>"
        #print tabl.style.background_gradient(cmap='summer')
        #sns.heatmap(tabl['Percentage of Interaction'], annot=True)
        Highest_value= tabl['Percentage of Interaction'][tabl['Percentage of Interaction']==tabl['Percentage of Interaction'].max()]
        Highest_value=Highest_value.to_dict()
        print "Highest percentage of interactions identified","<br/>"
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

All_Bonded_CSVfilename = '/tmp/'+ str(uuid.uuid4())+'_AllBondedCSV.csv'
Common_Bonded_CSVfilename = '/tmp/'+ str(uuid.uuid4())+'_CommonBondedCSV.csv'

print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>H - Bonded Interaction","</p>"
print '<p>Download: <a href=%s download>All Bonded</a>'% All_Bonded_CSVfilename
print "<a href=%s download>Common Bonded</a>" % Common_Bonded_CSVfilename ,"</p>"
print "<p align='center'>################################################################"  ,"</p>" 



#print Ligscompare

print "<button class='collapsible'>I. All Bonded Interactions - click here for basic statistical information </button>"#Start of click drop down
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"



    
print " Number of PDB IDs:", len(AllH_Lig_Resdict.keys()),"/", len(PDBID_LIST), "<br/>"

if bool(AllH_Lig_Resdict):
    print "[1: Interacting and 0: Non interacting]" , "<br/>"
    #print " Number of Ligand atoms:", Tot_Atom_count/len(finallistofAllatoms), "<br/>"
    percentage(AllH_Lig_Resdict) 

print "<br/>"    
print "Hydrogen bond distance"
if bool(AllH_distance_Lig_dict):
    distance_calc(AllH_distance_Lig_dict) 

print "</div>"#End of click drop down


print "<br/>"



print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""# initialization of first grid and column
if bool(AllH_Lig_Resdict):


    #print "All Hydrogen bondong details"
    bonded_csv= pd.DataFrame(AllH_Lig_Resdict)
    bonded_csv.to_csv(All_Bonded_CSVfilename)
    print pd.DataFrame.from_dict(AllH_Lig_Resdict).to_html(justify='center')

    #df=pd.DataFrame.from_dict(AllH_Lig_Resdict)
    #print df.style.applymap(lambda x: 'color: red' if pd.isnull(x) else '').to_html(justify='center')
else:
    print "All Hydrogen bondong details"
    print "<p No Hydrogen bonded interactions observed /p>"



print """
        </div>
    </div>
"""#closing of first column



###############Section for coloring amino acids based on physicochemical characeristics#################
#print All_graphdic_H
if bool(All_graphdic_H):
    
    templist4graph_H=[]
    
    graphdic1_All_H={}
    ###For all residues in selected PDB ids
    for k_all,v_all in All_graphdic_H.iteritems():
        for value_all in v_all:
            templist4graph_H.append(value_all)
            samp_all=sorted(list(set(templist4graph_H)))
        graphdic1_All_H.setdefault('%s'%k_all,[]).append(', '.join(samp_all))
        templist4graph_H=[]
    
    length_listofAllresidues_all_H=[]
    
    
    for key_all_H,value1_all_H in graphdic1_All_H.iteritems():    
        for i in value1_all_H:
            valu_all_H=i.split(', ')
            length_listofAllresidues_all_H.append(len(valu_all_H))
            
    length_ofcell_all_H=max(length_listofAllresidues_all_H)
    
    print """
    
       <div class="col-2-3">
         <div class="module">
       
    """# initialization of second column
    
    #print "<p style='font-size:20px; color:blue'> The All information of All Residues: H Bonded","</p>"
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'> List of all residues</th>"% length_ofcell_all_H
    print "</tr>"
    
    for key in sorted(graphdic1_All_H.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1_All_H in graphdic1_All_H[key]:
            dat1_All_H= g1_All_H.split(', ')
            for k3_H_All in dat1_All_H:
                print "<td align='center'>"
                #print k3
                if k3_H_All.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%k3_H_All
    
                if k3_H_All.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%k3_H_All                    
                
                if k3_H_All.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%k3_H_All
                    
                if k3_H_All.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%k3_H_All
                    
                if k3_H_All.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%k3_H_All
    
                if k3_H_All.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%k3_H_All
    
                if k3_H_All.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%k3_H_All
                
                print "</td>"
        print "<tr>"    
    print "</table>"

else:
    print "No All List for bonded interaction is available"
    
print """
        </div>
    </div>
"""#closing of first column

###################
#Script for generating web logo of All Ligand atoms 
#################
if bool(graphdic1_All_H):
    
    All_reviseddict_H={}
    All_reviseddict1_H={}
    for key in sorted(graphdic1_All_H):
        for i in graphdic1_All_H[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                All_reviseddict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in All_reviseddict_H.iteritems():
        counted=dict(Counter(n))
        All_reviseddict1_H.setdefault('%s'%m,{}).update(counted)
    
    ### generating WEBLOGO
    All_H_graph_filename = str(uuid.uuid4())
    All_H_graph_zipfilename='/tmp/'+All_H_graph_filename+'.zip'
    
    aminoacid_singlecode={}

    recoded={}
    for ligand_key, amino_frequency in All_reviseddict1_H.iteritems():
        #print ligand_key
        for i in ligand_key:
            
            for amino,frequency in amino_frequency.iteritems():
                
                for amino_3letter,code_frequency in aminoacid_code.iteritems():
                    if amino==amino_3letter:
                        recoded[code_frequency]=frequency
                        aminoacid_singlecode.setdefault('%s'%ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Frequency=1  
    instances=[]
    All_H_weblogo_collection=[]                  
    for ligand_key1, amino_frequency1 in aminoacid_singlecode.iteritems():
    
        for Amino1, number in amino_frequency1.iteritems():
    
            Frequency=1
            while Frequency <= number:
                instances.append(Seq(Amino1, IUPAC.protein))
                Frequency=Frequency+1
    
        m = motifs.create(instances)
    
        mymotif='/tmp/'+All_H_graph_filename+ligand_key1+'.svg'
        m.weblogo('%s'%mymotif,format='SVG',xaxis_label= '%s'%ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        All_H_weblogo_collection.append(mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in All_H_weblogo_collection)
    
    print """
       <div class="col-2-3">
         <div class="module">
    """# initialization of third column
    
    #print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    for All_H_image in sorted(All_H_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %All_H_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Nicot_image
        print "</div>"
    
    print "</div>"
    
    with ZipFile('%s'%All_H_graph_zipfilename, 'w') as All_H_myzip:
        for All_H_Images in All_H_weblogo_collection:
            
            All_H_myzip.write(All_H_Images)
else: 
    print "<p No Intercations!! /p>"
print """
        </div>
    </div>

</div>
"""#closing of first grid and third column

#print "<p align='center'>################################################################","</p>"
#print "<p style='font-size:20px; color:blue' align='center'>Common bonded interactions ","</p>"
#print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>II. Common bonded interactions - click here for basic statistical information </button>"
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"

print " Number of PDB IDs:", len(CommonH_Lig_Resdict.keys()),"/", len(PDBID_LIST), "<br/>"

if bool(CommonH_Lig_Resdict):
    print "Statistics of Bonded Intercations (1: Interacting and 0: Non interacting)" , "<br/>"
#print " Number of Ligand atoms:", Tot_Atom_count/len(finallistofAllatoms), "<br/>"
    percentage(CommonH_Lig_Resdict) 
else:
    print "<p No Common Bonded Intercations!! /p>"

print "<br/>"
print "Hydrogen bond distance"
if bool(CommonH_distance_Lig_dict):
    distance_calc(CommonH_distance_Lig_dict)

else:
    print "<p No Common Bonded Intercations!! /p>"
    

print "</div>" #close of collapsable section


print "<br/>"
print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""#initializing second grid and first column



########For Common amino residues
if bool(CommonH_Lig_Resdict):

    Com_bonded_csv = pd.DataFrame(CommonH_Lig_Resdict)
    Com_bonded_csv.to_csv(Common_Bonded_CSVfilename)
   # print "Common Hydrogen bondong details"

    print pd.DataFrame.from_dict(CommonH_Lig_Resdict).to_html(justify='center')

else:
    print "<p No Common Bonded Intercations!! /p>"
print """
        </div>
    </div>


"""#closing of first column in second grid
templist4graph_H=[]

graphdic1_Common_H={}
if bool(common_graphdic_H):
    

    ###For all residues in selected PDB ids
    for k_all,v_all in common_graphdic_H.iteritems():
        for value_all in v_all:
            templist4graph_H.append(value_all)
            samp_all=sorted(list(set(templist4graph_H)))
        graphdic1_Common_H.setdefault('%s'%k_all,[]).append(', '.join(samp_all))
        templist4graph_H=[]
    
    length_listofAllresidues_common_H=[]
    
    
    for key_common_H,value1_common_H in graphdic1_Common_H.iteritems():    
        for i in value1_common_H:
            valu_common_H=i.split(', ')
            length_listofAllresidues_common_H.append(len(valu_common_H))
            
    length_ofcell_common_H=max(length_listofAllresidues_common_H)
    print """
    
       <div class="col-2-3">
         <div class="module">
       
    """# initialization of secon grid second column
    
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>List of common residues</th>"% length_ofcell_common_H
    print "</tr>"
    for key in sorted(graphdic1_Common_H.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1_Common_H in graphdic1_Common_H[key]:
            dat1_Common_H= g1_Common_H.split(', ')
            for k3_H_Common in dat1_Common_H:
                print "<td align='center'>"
                #print k3
                if k3_H_Common.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%k3_H_Common
    
                if k3_H_Common.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%k3_H_Common                    
                
                if k3_H_Common.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%k3_H_Common
                    
                if k3_H_Common.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%k3_H_Common
                    
                if k3_H_Common.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%k3_H_Common
    
                if k3_H_Common.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%k3_H_Common
    
                if k3_H_Common.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%k3_H_Common
                
                print "</td>"
        print "<tr>"    
    print "</table>"

else:
    print "<p No All Common Bonded Intercations for these structures!! /p>"
    
print """
        </div>
    </div>
"""#closing of second column

#web log for common hydrogen
if bool(graphdic1_Common_H):
    
    Common_reviseddict_H={}
    Common_reviseddict1_H={}
    for key in sorted(graphdic1_Common_H):
        for i in graphdic1_All_H[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Common_reviseddict_H.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in Common_reviseddict_H.iteritems():
        counted=dict(Counter(n))
        Common_reviseddict1_H.setdefault('%s'%m,{}).update(counted)
    
    ### generating WEBLOGO
    Common_H_graph_filename = str(uuid.uuid4())
    Common_H_graph_zipfilename='/tmp/'+Common_H_graph_filename+'.zip'
    
    aminoacid_singlecode={}
   # aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    #     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     #    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
      #   'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    
    for ligand_key, amino_frequency in Common_reviseddict1_H.iteritems():
        #print ligand_key
        for i in ligand_key:
            
            for amino,frequency in amino_frequency.iteritems():
                
                for amino_3letter,code_frequency in aminoacid_code.iteritems():
                    if amino==amino_3letter:
                        recoded[code_frequency]=frequency
                        aminoacid_singlecode.setdefault('%s'%ligand_key,{}).update(recoded)
                        recoded={}
    
    
    Frequency=1  
    instances=[]
    Common_H_weblogo_collection=[]                  
    for ligand_key1, amino_frequency1 in aminoacid_singlecode.iteritems():
    
        for Amino1, number in amino_frequency1.iteritems():
            Frequency=1
            while Frequency <= number:
                instances.append(Seq(Amino1, IUPAC.protein))
                Frequency=Frequency+1
    
        m = motifs.create(instances)
    
        mymotif='/tmp/'+All_H_graph_filename+ligand_key1+'.svg'
        m.weblogo('%s'%mymotif,format='SVG',xaxis_label= '%s'%ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Common_H_weblogo_collection.append(mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in Common_H_weblogo_collection)
    
    print """
       <div class="col-2-3">
         <div class="module">
    """# initialization of third column
    
    #print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    for All_H_image in sorted(Common_H_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %All_H_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Nicot_image
        print "</div>"
    
    print "</div>"
    
    with ZipFile('%s'%Common_H_graph_zipfilename, 'w') as Common_H_myzip:
        for Common_H_Images in All_H_weblogo_collection:
            
            Common_H_myzip.write(Common_H_Images)
else:
    print "<p No Intercations!! /p>"

print """
        </div>
    </div>

</div>
"""#closing of second grid and third column

All_NonBonded_CSVfilename = '/tmp/'+ str(uuid.uuid4())+'_AllNonBondedCSV.csv'
Common_NonBonded_CSVfilename = '/tmp/'+ str(uuid.uuid4())+'_CommonNonBondedCSV.csv'

print "<p align='center'>################################################################","</p>"
print "<p style='font-size:20px; color:blue' align='center'>Non-bonded Interaction","</p>"
print "<p>Download: <a href=%s download>All Non-bonded</a>" % All_NonBonded_CSVfilename
print "<a href=%s download>Common Non-bonded</a>" % Common_NonBonded_CSVfilename ,"</p>" 
print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>I. All non-bonded interactions - click here for basic statistical information </button>"
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"


print " Number of PDB IDs:", len(AllNH_Lig_Resdict.keys()),"/", len(PDBID_LIST), "<br/>"

if bool(AllNH_Lig_Resdict):
    print "Statistics of Bonded Intercations (1: Interacting and 0: Non interacting)" , "<br/>"
    #print " Number of Ligand atoms:", Tot_Atom_count/len(finallistofAllatoms), "<br/>"
    percentage(AllNH_Lig_Resdict)

print "<br/>"
print " Contact distance"
if bool(AllNH_distance_Lig_dict):
    distance_calc(AllNH_distance_Lig_dict)


print "</div>"#close of collapsable section


print "<br/>"
print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""#initializing third grid and first column
if bool(AllNH_Lig_Resdict):
    #print "All Non-bondong details"
    Non_Allbonded_csv = pd.DataFrame(AllNH_Lig_Resdict)
    Non_Allbonded_csv.to_csv(All_NonBonded_CSVfilename)

    print pd.DataFrame.from_dict (AllNH_Lig_Resdict).to_html(justify='center') 
else: 
    print "<p No Non-bonded Interaction in these structures"

print """
        </div>
    </div>


"""#closing of first column in third grid
templist4graph_NH=[]

graphdic1_All_NonH={}
if bool(All_graphdic_NonH):
    

    ###For all residues in selected PDB ids
    for k_all,v_all in All_graphdic_NonH.iteritems():
        for value_all in v_all:
            templist4graph_NH.append(value_all)
            samp_all=sorted(list(set(templist4graph_NH)))
        graphdic1_All_NonH.setdefault('%s'%k_all,[]).append(', '.join(samp_all))
        templist4graph_NH=[]
    
    length_listofAllresidues_all_NonH=[]
    
    
    for key_all_H,value1_all_H in graphdic1_All_NonH.iteritems():    
        for i in value1_all_H:
            valu_all_H=i.split(', ')
            length_listofAllresidues_all_NonH.append(len(valu_all_H))
            
    length_ofcell_all_NonH=max(length_listofAllresidues_all_NonH)
    
    print """
    
       <div class="col-2-3">
         <div class="module">
       
    """# initialization of second column in third grid
    
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>List of all residues</th>"% length_ofcell_all_NonH
    print "</tr>"
    
    for key in sorted(graphdic1_All_NonH.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1_All_NonH in graphdic1_All_NonH[key]:
            dat1_All_NonH= g1_All_NonH.split(', ')
            for k3_NonH_All in dat1_All_NonH:
                print "<td align='center'>"
                if k3_NonH_All.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%k3_NonH_All
    
                if k3_NonH_All.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%k3_NonH_All                    
                
                if k3_NonH_All.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%k3_NonH_All
                    
                if k3_NonH_All.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%k3_NonH_All
                    
                if k3_NonH_All.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%k3_NonH_All
    
                if k3_NonH_All.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%k3_NonH_All
    
                if k3_NonH_All.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%k3_NonH_All
                
                print "</td>"
        print "<tr>"    
    print "</table>"

else:
    print "<p No All list of nonbonded residues /p>"
    
    
print """
        </div>
    </div>
"""#closing of first column

###################
#Script for generating web logo of All Ligand atoms 
#################
if bool(graphdic1_All_NonH):
    
    All_reviseddict_NonH={}
    All_reviseddict1_NonH={}
    for key in sorted(graphdic1_All_NonH):
        for i in graphdic1_All_NonH[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                All_reviseddict_NonH.setdefault('%s'%key,[]).append(se[0])
    
    for m,n in All_reviseddict_NonH.iteritems():
        counted=dict(Counter(n))
        All_reviseddict1_NonH.setdefault('%s'%m,{}).update(counted)
        
    
    ### generating WEBLOGO
    All_NonH_graph_filename = str(uuid.uuid4())
    All_NonH_graph_zipfilename='/tmp/'+All_NonH_graph_filename+'.zip'
    
    aminoacid_singlecode={}
    #aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
      #   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
       ## 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    
    for ligand_key, amino_frequency in All_reviseddict1_NonH.iteritems():
        #print ligand_key
        for i in ligand_key:
            
            for amino,frequency in amino_frequency.iteritems():
                
                for amino_3letter,code_frequency in aminoacid_code.iteritems():
                    if amino==amino_3letter:
                        recoded[code_frequency]=frequency
                        aminoacid_singlecode.setdefault('%s'%ligand_key,{}).update(recoded)
                        recoded={}
    
    Frequency=1  
    instances=[]
    All_NonH_weblogo_collection=[]                  
    for ligand_key1, amino_frequency1 in aminoacid_singlecode.iteritems():
    
        for Amino1, number in amino_frequency1.iteritems():
    
            Frequency=1
            while Frequency <= number:
                instances.append(Seq(Amino1, IUPAC.protein))
                Frequency=Frequency+1
    
        m = motifs.create(instances)
        
        mymotif='/tmp/'+All_NonH_graph_filename+ligand_key1+'.svg'
        m.weblogo('%s'%mymotif,format='SVG',xaxis_label= '%s'%ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        All_NonH_weblogo_collection.append(mymotif)
        instances=[]
    weblogo_images=' '.join(str(x) for x in All_NonH_weblogo_collection)
    
    print """
       <div class="col-2-3">
         <div class="module">
    """# initialization of third column
    
    
    #print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    for All_NonH_image in sorted(All_NonH_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %All_NonH_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Nicot_image
        print "</div>"
    
    print "</div>"
    
    with ZipFile('%s'%All_NonH_graph_zipfilename, 'w') as All_NonH_myzip:
        for All_NonH_Images in All_NonH_weblogo_collection:
            
            All_NonH_myzip.write(All_NonH_Images)

else:
    print "<p No Intercations!! /p>"
print """
        </div>
    </div>

</div>
"""#closing of third grid and third column
    

#print "<p align='center'>################################################################","</p>"
#print "<p style='font-size:20px; color:blue' align='center'>Common Nonbonding Interaction","</p>"
#print "<p align='center'>################################################################"  ,"</p>"

print "<button class='collapsible'>II. Common non-bonded interactions - click here for basic statistical information </button>"
print "<div class='contentsection'>"

print "<p style='font-size:20px; color:black' align='center'>"


print " Number of PDB IDs:", len(CommonNH_Lig_Resdict.keys()),"/", len(PDBID_LIST), "<br/>"

if bool(CommonNH_Lig_Resdict):
    print "Statistics of Bonded Intercations (1: Interacting and 0: Non interacting)" , "<br/>"
    #print " Number of Ligand atoms:", Tot_Atom_count/len(finallistofAllatoms), "<br/>"
    percentage(CommonNH_Lig_Resdict)

print "<br/>"
print "Contact distance"
if bool(CommonNH_distance_Lig_dict):
    distance_calc(CommonNH_distance_Lig_dict)

print "</div>"#close of collapsable section


print "<br/>"
print """
<div class="grid">
   <div class="col-2-3">
     <div class="module">
   
"""#initializing fourth grid and first column
if bool(CommonNH_Lig_Resdict):
    #print "Common Non-bondong details"
    Common_Nonbonded_csv = pd.DataFrame(CommonNH_Lig_Resdict)
    Common_Nonbonded_csv.to_csv(Common_NonBonded_CSVfilename)
    print pd.DataFrame.from_dict(CommonNH_Lig_Resdict).to_html(justify='center')

else:
    print "<p No Common non-bonded interactions! /p>"


print """
        </div>
    </div>


"""#closing of first column in third grid
##### Colored amino acid table#########
templist4graph_NH=[]

graphdic1_Common_NonH={}
if bool(common_graphdic_NonH):
    

    ###For all residues in selected PDB ids
    for k_all,v_all in common_graphdic_NonH.iteritems():
        for value_all in v_all:
            templist4graph_NH.append(value_all)
            samp_all=sorted(list(set(templist4graph_NH)))
        graphdic1_Common_NonH.setdefault('%s'%k_all,[]).append(', '.join(samp_all))
        templist4graph_NH=[]
    
    length_listofAllresidues_Common_NonH=[]
    
    
    for key_all_H,value1_all_H in graphdic1_Common_NonH.iteritems():    
        for i in value1_all_H:
            valu_all_H=i.split(', ')
            length_listofAllresidues_Common_NonH.append(len(valu_all_H))
            
    length_ofcell_Common_NonH=max(length_listofAllresidues_Common_NonH)
    
    print """
    
       <div class="col-2-3">
         <div class="module">
       
    """# initialization of second column in third grid
    
    print "<table border='1'>"
    print "<tr>"
    print "<th col width='60'>Ligand Atoms</th>" 
    print "<th  colspan='%d'>List of all residues</th>"% length_ofcell_Common_NonH
    print "</tr>"
    
    for key in sorted(graphdic1_Common_NonH.iterkeys()):
        print "<td align='center'>%s</td>" %key
        for g1_Common_NonH in graphdic1_Common_NonH[key]:
            dat1_Common_NonH= g1_Common_NonH.split(', ')
            for k3_NonH_Common in dat1_Common_NonH:
                print "<td align='center'>"
                if k3_NonH_Common.startswith(('ALA','ILE','LEU','MET','MSE','VAL')):                
                    print "<b><font color='pink'>%s</font></b>"%k3_NonH_Common
    
                if k3_NonH_Common.startswith(('PHE','TRP', 'TYR')):                
                    print " <b><font color='orange'>%s</font></b>"%k3_NonH_Common                    
                
                if k3_NonH_Common.startswith(('LYS','ARG', 'HIS')):                
                    print " <b><font color='red'>%s</font></b>"%k3_NonH_Common
                    
                if k3_NonH_Common.startswith(('GLU','ASP')):                
                    print " <b><font color='green'>%s</font></b>"%k3_NonH_Common
                    
                if k3_NonH_Common.startswith(('ASN','GLN','SER','THR')):                
                    print " <b><font color='blue'>%s</font></b>"%k3_NonH_Common
    
                if k3_NonH_Common.startswith(('GLY','PRO')):                
                    print " <b><font color='magenta'>%s</font></b>"%k3_NonH_Common
    
                if k3_NonH_Common.startswith(('CYS','CME')):                
                    print " <b><font color='yellow'>%s</font></b>"%k3_NonH_Common
                
                print "</td>"
        print "<tr>"    
    print "</table>"

else:
    print "<p No Compiles list of Non-bonded Interactions /p>"
print """
        </div>
    </div>
"""#closing of first column

###################
#Script for generating web logo of All Ligand atoms 
#################
Common_reviseddict_NonH={}
Common_reviseddict1_NonH={}
if bool(graphdic1_Common_NonH):

    for key in sorted(graphdic1_Common_NonH):
        for i in graphdic1_Common_NonH[key]:
            tems=i.split(', ')
            for items in tems:
                se=re.split('([0-9])' , items)
                Common_reviseddict_NonH.setdefault('%s'%key,[]).append(se[0])
    
    
    for m,n in Common_reviseddict_NonH.iteritems():
        counted=dict(Counter(n))
        Common_reviseddict1_NonH.setdefault('%s'%m,{}).update(counted)
        
    ### generating WEBLOGO
    Common_NonH_graph_filename = str(uuid.uuid4())
    Common_NonH_graph_zipfilename='/tmp/'+Common_NonH_graph_filename+'.zip'
    
    aminoacid_singlecode={}
    #aminoacid_code={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     #    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
      #   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
       #  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    recoded={}
    
    for ligand_key, amino_frequency in Common_reviseddict1_NonH.iteritems():
    
        for i in ligand_key:
            for amino,frequency in amino_frequency.iteritems():
                
                for amino_3letter,code_frequency in aminoacid_code.iteritems():
                    if amino==amino_3letter:
                        recoded[code_frequency]=frequency
                        aminoacid_singlecode.setdefault('%s'%ligand_key,{}).update(recoded)
                        recoded={}
    
    Frequency=1  
    instances=[]
    Common_NonH_weblogo_collection=[]                  
    for ligand_key1, amino_frequency1 in aminoacid_singlecode.iteritems():
        for Amino1, number in amino_frequency1.iteritems():
            Frequency=1
            while Frequency <= number:
                instances.append(Seq(Amino1, IUPAC.protein))
                Frequency=Frequency+1
    
        m = motifs.create(instances)
    
    
        mymotif='/tmp/'+Common_NonH_graph_filename+ligand_key1+'.svg'
        m.weblogo('%s'%mymotif,format='SVG',xaxis_label= '%s'%ligand_key1,show_errorbars= False, color_scheme= 'color_chemistry')
        Common_NonH_weblogo_collection.append(mymotif)
        instances=[]
    
    weblogo_images=' '.join(str(x) for x in Common_NonH_weblogo_collection)
    
    print """
       <div class="col-2-3">
         <div class="module">
    """# initialization of third column
    
    #print "<p style='font-size:20px; color:brown'> Weblogo showing the frequency of Residues binding to Ligand atoms for the selected structures:"
    
    print "<div class='weblogo_row'>"
    
    for Common_NonH_image in sorted(Common_NonH_weblogo_collection):
        print "<div class='weblogo_column'>"
        print "<embed src='%s#page=1&view=FitH ' />" %Common_NonH_image
        #print "<iframe src='%s#page=1&view=FitH ' width='200' height='100' border='0'></iframe>"%Nicot_image
        print "</div>"
    
    print "</div>"
    
    
    with ZipFile('%s'%Common_NonH_graph_zipfilename, 'w') as Common_NonH_myzip:
        for Common_NonH_Images in Common_NonH_weblogo_collection:
            
            Common_NonH_myzip.write(Common_NonH_Images)
else:
    print "<p No Intercations! /p>"
print """
        </div>
    </div>

</div>
"""#closing of fourth grid and fourth column
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
