#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 10:30:10 2018

@author: xhasam
"""

import cgi, cgitb
import urllib
import urllib2
import requests
from collections import Counter
import re
print 'Content-type:text/html\r\n\r\n'
print '<html>'
print '<head>'
print """
<script type='text/javascript' src='jsmol/jsmol/JSmol.min.js'></script>

	<script type='text/javascript'>
	var Info = {
	  width: 200,
	  height: 100,
	  j2sPath: 'jsmol/jsmol/j2s',
	}
	</script>
"""
atomssele=['C1B', 'C1D', 'C2A', 'C2D', 'C2N', 'C3B', 'C3D', 'C3N', 'C4D', 'C4N', 'C5A', 'C5B', 'C5D', 'C5N', 'C6A', 'C6N', 'C7N', 'C8A', 'N1A', 'N1N', 'N3A', 'N6A', 'N7A', 'N7N', 'O1A', 'O1N', 'O1X', 'O2A', 'O2B', 'O2D', 'O2N', 'O2X', 'O3', 'O3B', 'O3D', 'O3X', 'O4B', 'O4D', 'O5B', 'O7N', 'P2B', 'PN']
#atomssele=['C2O','O1D','O2A']
atomseen=set()
resseen=set()
linkspdb=['https://files.rcsb.org/view/3wxb.pdb','https://files.rcsb.org/view/3p19.pdb']
residuelist=['GLY13A','ILE157A', 'LEU66A', 'ARG41A', 'VAL68A','MET215A', 'TYR178A','ILE18A','ASN15A', 'ARG16A','ASN96A']
#links=['http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=3wxb&ligtype=01&ligno=01','http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=3p19&ligtype=01&ligno=01']
links=['http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=3wxb&ligtype=01&ligno=01','http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=3p19&ligtype=01&ligno=01']
pdbids=['3wxb','3p19']
dict1={}
dict2={}
listsel_resname=[]
listresnum=[]
listreschain=[]
listfinal_res_to_extract=[]
lig_detail_to_extract=[]
lineofinterest=[]
the_list=[]
#for res in residuelist:
 #   print res
  #  print re.split('([0-9])' , res)

for link,pdblink,pdbid in zip(links,linkspdb,pdbids):

    #for pdblink in linkspdb:
    #print link
    weblink=requests.get(link, stream=True)
    #print weblink
    for atomlines in weblink.iter_lines():
        atomlines1=atomlines.strip()
        #print atomlines1
        if atomlines1.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
            sel=atomlines1.split()
            res=sel[3]
            resno=sel[4]
            resnum=res+resno
            ligatm=sel[8]
            ligname=sel[9]
            lignumber=sel[10]
            ligchain=sel[11]
            ligdetail=ligname + '_' + lignumber + '_' + ligchain
            ligdetail1=ligname + '_' + ligchain+lignumber
            reschain=sel[5]
            #print reschain
            #print ligatm,resnum
            if ligatm in atomssele:
                sel_resname=res
                sel_resnum=resno
                sel_reschain=reschain
                final_res_to_extract=sel_resname + '_' + sel_resnum + '_' + sel_reschain
                #print final_res_to_extract
                listfinal_res_to_extract.append(final_res_to_extract)
                lig_detail_to_extract.append(ligdetail)
                #print listfinal_res_to_extract
                set_listfinal=list(set(listfinal_res_to_extract))
                #set_lig_final=list(set(lig_detail_to_extract))
            
    #print "LIGDETAILS", ligdetail
    #print "LIGDETAILS", ligdetail1
    ligsplit=ligdetail.split('_')
    ligsplit1=ligdetail1.split('_')
    print ligsplit[1] 
    ligcode=ligsplit[0]+ligsplit[2]+ligsplit[1]
    filecode=pdbid +ligcode+'.pdb'
    #print ligsplit1              
    print "SET RESDUES",set_listfinal
    #print len(set_listfinal)
    the_list.append(set_listfinal)

    #print linkspdb
    #print "LIGAND",set_lig_final
    #for pdblink in linkspdb:
    print "PDB LINK",pdblink
    print "PDB LINK",link
    with open(filecode,'w') as fout:
        weblink1=requests.get(pdblink, stream=True)
        for residuelines in weblink1.iter_lines():
            if residuelines.startswith("ATOM"):
                lines1=residuelines.split()
                #print lines1
                for residue in set_listfinal:
                    
                    
                    ressplit=residue.split('_')
                    if ressplit[0]==lines1[3] and ressplit[2] == lines1[4] and ressplit[1] == lines1[5] :
                        reslines= '   '.join(lines1)
                        
                        print reslines
                        
                        fout.write(reslines+'\n')
            if residuelines.startswith("HETATM"):
                lines1=residuelines.split()
                ligsplit=ligdetail.split('_')
                if ligsplit[0]==lines1[3] and ligsplit[2] == lines1[4] and ligsplit[1] == lines1[5] :
                    liglines= '   '.join(lines1)
                    print liglines
                    fout.write(liglines+'\n')
            if residuelines.startswith("HETATM"):
                lines1=residuelines.split()
                ligsplit1=ligdetail1.split('_')
                if ligsplit1[0]==lines1[3] and ligsplit1[1] == lines1[4] :
                    liglines1= '   '.join(lines1)
                    print liglines1
                    fout.write(liglines1+'\n')
    print filecode
    print """
    <script type='text/javascript'>
    jmolApplet0 = Jmol.getApplet('jmolApplet0', Info);
    jmolRadio(jmolApplet0,'background black; load https://files.rcsb.org/view/%s; wireframe only')</script>
    """% filecode        
    #print "check",set_listfinal                  #checkline=lines1
    listfinal_res_to_extract=[]
    set_listfinal=[]

listpdb=['1bdm','3aek']

#print the_list
#for types in the_list:
#    print types
##    print len(types)
#    for names in types:
#        print names
    
                                #if ressplit[0] and ressplit[1] and ressplit[2] in residuelines:
                                    #print residuelines
                        #if  sel_resname and sel_resnum and sel_reschain in lines1:
                         #   print line

print "</body>"
print "</html>"                    
                
                
