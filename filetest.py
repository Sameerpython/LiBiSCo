#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 14:38:31 2018

@author: xhasam
"""

#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:25:22 2018

@author: xhasam
"""
from itertools import chain
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

pd.set_option('display.max_colwidth', -1)



print "Content-type:text/html\r\n\r\n"
print '<html>'
print '<head>'
print '<title>Hello Word - First CGI Program</title>'
print '</head>'
print '<body>'

########################
atom_C=sorted(['ND','C1D','C2D','C3D','C4D','CMD'])
mydictcheck={'3aek': ['http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=3aek&ligtype=02&ligno=01'], '2ynm': ['http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=2ynm&ligtype=06&ligno=01']}

lresidue=[]
latom=[]
ldistance=[]
residue={}
atmname={}
dicresidue_unique={}
residue_seen=set()
atom_seen=set()
finalset=set()
combines_listdata=[]
ATMA_listdata=[]

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
                    if line.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                            line=line.split()
                            ligname=line[9]
                            atm=line[8]
                            res=line[3]
                            residuenum=line[4]
                            distance=line[12]
                            resnum=res+residuenum
#                               #appending each residue and its position to list called lresidue
                            lresidue.append(resnum)
                            #appending each ligand atom to list called latom
                            latom.append(atm)
                            #appending distance of each interaction to ldistance
                            ldistance.append(distance)
                            #creating a set for residue with position
                            residue_seen.add(resnum)
                            #creating a set for each ligand atom
                            atom_seen.add(atm)
                            #making a dictionary with list comtaining residue name and position
                            residue.setdefault('%s'%id,[]).append(resnum)
                            #making a dictionary with list comataing ligand atoms
                            atmname.setdefault('%s'%id,[]).append(atm)

#converting set to a list for residue_seen
   list_residue_seen=list(residue_seen)
#converting set to alist for atom_seen
   list_atom_seen=list(atom_seen)
#merging the two lists:lresidue and latom based on their index position
   merged_atm_residue=[x1 +"        "+ y1 for x1, y1 in zip(latom,lresidue)]

#merging the three lists:lresidue, latom and ldistance based on their index position
   merged_atm_residue_distance=[x2 +"&nbsp "+ "&nbsp"+ y2 +"&nbsp "+ "&nbsp"+ z2 for x2, y2, z2 in zip(latom,lresidue,ldistance)]

#sorting the merged list based on aton name
   sort_merged_atm_residue=sorted(merged_atm_residue)

#removing duplicated from sort_merged_atm_residue
   for word in sort_merged_atm_residue:
           if word not in finalset:
                   finalset.add(word)
#converting finalset to a list function
   list_finalset=list(finalset)
   lig_tabledic={}
   appended_lig_tabledic={}
   #print "RING C", id,lig 

   ATOM=[i for e in atom_C  for i in list_finalset if e in i]   
   for PC in ATOM:
       subC=PC.split()
       ATMA_listdata.extend(subC)
   print ATMA_listdata
   size=int(len(ATMA_listdata)/2)
   df_F1=pd.DataFrame(np.array(ATMA_listdata).reshape(size,2),columns=['LIGAND ATOM','RESIDUES FOR PDB %s'%id])
   groupedF1=df_F1.groupby(['LIGAND ATOM'])['RESIDUES FOR PDB %s'%id].apply(lambda x: ', '.join(x.astype(str))).reset_index()
   print "<td>" 
   print groupedF1.to_html(justify='center') 
   lig_dict=groupedF1.to_dict('split')
   lig_dict_new={}
   for ligs in lig_dict.iterkeys():
       if ligs=='data':
           for atms in lig_dict[ligs]:
               key, value=atms[0], atms[1:]
               lig_dict_new[key]=value
   #print lig_dict_new
   
               #lig_tabledic={'%s'%id:lig_dict_new}
               #print lig_tabledic
   #for dickey, dicvalue in appended_lig_tabledic:
               appended_lig_tabledic['%s'%id]=lig_dict_new 
   print "check update", appended_lig_tabledic     
   #appended_lig_tabledic.update(lig_tabledic) 
   #print "check update", appended_lig_tabledic            
               
   print "</td>"

   lresidue=[]
   latom=[]
   ldistance=[]
   combines_listdata=[]
   ATMA_listdata=[]
   residue_seen.clear()
   finalset.clear()
   








print '</body>'
print '</html>'


            
        