#!/usr/bin/python2.7

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
if form.getvalue('dropdown'):
   text_content = form.getvalue('dropdown')
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

print "<title>CoFactComp</title>"
print "</head>"
print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
print "<div align='center'>"
print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
print "</div>"
print "<body>"
print "<ul>"
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href=''>Contact</a></li>"
print "</ul>"
print "<div align='center'>"
print "<h2> The Selected Cofactor for comparison is : %s</h2>" % text_content
#print "<p> Select ligand of interest for the PDB IDs listed to compare the binding sites</p>"

#capturing the entered pdb ids into list
#print text_content

print "<p> Enter PDB ids bound to %s here seprated by comma for comparison (For eg:3WXB,3P19):</p>"% text_content
print  "<form action='%s.py' method = 'post' target = '_blank'>" % text_content
print  "<textarea rows='4' cols='50' name = 'textcontent' cols = '40' rows = '40'>"
print  "</textarea>"
print  "<p><input type = 'submit' value = 'Submit' /></p>"
print  "</form>" 


print "<div class='footer'>"
print "<p><img src='University.jpg' float='left'  width='40' height='40'<a href='https://bioenv.gu.se/english/staff?userId=xarohe'>Prof. Henrik Aronsson Group</a>, Department of Biological and Environmental Sciences, University of Gothenburg, Sweden. </p>"

print "</div>"
print "</div>"




print "</body>"
print "</html>"
