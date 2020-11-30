#!/usr/bin/python

# Import modules for CGI handling 
import cgi, cgitb 
import urllib
import requests
import re
import sys, os
url="https://files.rcsb.org/view/"


# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from field
if form.getvalue('textcontent'):
   text_content = form.getvalue('textcontent')
else:
   text_content = "Not entered"

print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"

#HTML style details for web page
print "<style>"
print "ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}"
print "li{float:left;}"
print "li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}"
print "li a:hover { background-color: #111111;}"
print "table, th, td { border: 2px solid black;}"
print ".footer { position: absolute; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }"
print "</style>"
#Style ends here

print "<title>LiBiSCo</title>"
print "</head>"
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
print "<div align='center'>"
print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
print "</div>"
print "<body>"
print "<ul>"
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href=''>Contact</a></li>"
print "</ul>"
print "<div align='center'>"
print "<h2> Below is shown Bound Ligand information for the PDB IDs entered: %s</h2>" % text_content
print "<p> Select ligand of interest for the PDB IDs listed to compare the binding sites</p>"

#Capturing URL address to get the ligand name for further processing
url1 = os.environ['HTTP_HOST']
uri = os.environ['REQUEST_URI']
filename= url1+uri
take=filename.split('/')
take2= take[-1]
take3=take2.split('.')
final_filename= take3[0]
#URL address capturing finishes

#capturing the entered pdb ids into list
for i in text_content:
        text_content1=text_content.replace(' ', '')
        l=text_content1.split(',')
f2_list=[]
combined_list=[]
#print "entered id's",l
pdbid_list=[]
url_list=[]
mydict={}
count=1
for i in l:
#        print "looping over the entered pdb ids:",i
        each= i # i am removing the +"="+"[]"
        pdbid_list.append(each)
#print "list generated for individual PDB id entered:", pdbid_list

#linking the PDB url address to the selected PDB ids
for x in l:
#        print "linking each id to its url address(below each id is its url address):",x
#        print "<br/>"
        pdbid=url+x+".pdb"
#        print pdbid
        url_list.append(pdbid)
#print "url list content is :",url_list

#Listing the PDB ids and their bound ligands
for link,id  in zip(url_list,pdbid_list):
#       print "This is the %d st/nd number of url %s"% (count,link)
        count=count+1
        f=urllib.urlopen(link)
        f=f.readlines()
#       print f
#!      f='\n'.join(f)
        for e in f:
#               print "this is:",e.split()
                if e.startswith('FORMUL'):
                        e= e.split()
        #               print e
                        if e[2]!="HOH":
                                f2=e[2]
                        #       print "ligands extacted %s from %s:"% (f2,link)
#                               id_list.append(f2)
#                               print id_list
                                mydict.setdefault('%s'%id,[]).append(f2)
#print "<br/>"
#print "Dictionary od PDBIDS", mydict

# Form for selecting the ligands for interest for the USER
print "<form action='%s_Result.py' method = 'post' target = '_blank'>"% final_filename 
print "<table style=width:50%>"
print "<tr>"
print "<th>PDB ID</th>"
print "<th colspan=2>LIGANDS</th>"
for k,dk in mydict.iteritems():
	count= len(dk)
	print "<tr>"
        print "<th rowspan='%d'>"% count,k,"</th>"
        #print "</tr>"
        for x in dk:
                print "<td align=center>", "%s" % x,"</td>"
                print "<td align=center>"
		print  "<input type='checkbox' name='%s' value='%s'/>" % (k,x)
		print "</td>"
		print "</tr>"
#                print "<br/>"
#        print "<br/>"
print "</table>"
print "<br/>"
print "<input type = 'submit'  value = 'Submit'  />"
print "</form>"


 



print "</div>"




print "</body>"
print "</html>"
