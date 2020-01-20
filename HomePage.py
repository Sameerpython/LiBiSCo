#!/usr/local/Anaconda2.7/bin/python2.7

# Import modules for CGI handling 
import cgi, cgitb 
# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from fields
if form.getvalue('subject'):
   subject = form.getvalue('subject')
else:
   subject = "Not set"


print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"

print "<style>"
#header styling
print"""
.images{
    width: 100%;
    height: 500px;
    background-image: url('/opt/lampp/htdocs/LiBiSCo/Pic.gif') ;
    background-size: cover;
    padding-left: 80px;
    
}
}
"""
#########

print "div.container { width:100%; border: 1px solid grey;}"
print "ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}"
print "li{float:left;}"
print "li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}"
print "li a:hover { background-color: #111111;}"
print ".footer { position: fixed; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }"
print " #container1{ width:1000px; height=1500px; line-height:1.6;}"

#style for divindg into 2 columns

print "* {box-sizing: border-box;}"
print ".column {float: left;width: 50%;padding: 10px;height: 300px;}"
print ".row:after {content: "";display: table;clear: both;}"


print "</style>"

print "<title>LiBiSCo</title>"
print "</head>"
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
#print "<div align='center'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "</div>"

#print "<div id='container3' style='position:relative;'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "<h1 style='position:absolute; top:100px; left:20px;'>LiBisCo </h1>  "
#print "<h1 style='position:absolute; top:100px; left:20px;'> Ligand Binding Site Comparison </h1>   " 
#print "</div>"

print "<div class='images'>"
#print "<div class='main'>"
#print "<style='font-size:120px'>LiBiSCo</font>","</br>"
#print "<font color='red'> Li </font>gand <font color='red'>Bi</font>nding <font color='red'>S</font>ite <font color='red'>Co</font>mparison"
print"  </div>"
#print "</div>"



print "<body>"
print "<ul>"
print "<li><a href=''>Home</a></li>"
print "<li><a href='Contact.py'>Contact</a></li>"
print "</ul>"
#Body size Determination


print "<div align='center'>"
print "<h2> Comparing the Binding Residues for the Selected Protein Structures</h2>"
#print "</br>"
print "<div id='container1'>"
print "<p align='justify'>Proteins bind to ligands to perform all kinds of important and essential cellular processes. They bind via a network of weak, noncovalent intermolecular interactions such as hydrogen bonding, hydrophobic and electrostatic interactions. Thus, binding of substrate is required for many proteins to function properly. Ligands are recognized in the binding pockets of proteins through surface exposed properties of amino acids. Ligands conformational flexibility make them bind to homologous proteins as well to proteins with remote similarity with reference to the cavity shape and folding pattern. These proteins despite of their divergent shapes in active site are recognized by superposition based on their similarity in key pharmacophoric features.</p>"
print "<p align='justify'> However, there are proteins that are different in fold as well as cavity shape which are not recognized or detected by these programs. In such scenarios, it is important to compare residues that interact with these structurally similar ligands among proteins that are divergent in their fold as well in cavity shape. The program <b> CoFactComp </b> helps in comparing the binding site residues among proteins that bind to structurally similar ligands irrespective of having structural similarity or not at the fold level or cavity shape.</p>" 

print "<div class='row'>"
print "<div class='column'>"

print "<h2> Ligand Based Analysis</h1>"

print "<p> Enter PDB ids here separated by comma (eg:3WXB,3P19):</p>"
print  "<form action='LigPage.py' method = 'post' target = '_blank'>"
print  "<textarea rows='4' cols='50' name = 'textcontent' cols = '40' rows = '40'>"
print  "</textarea>"
#print  "<p>Enter Ligand Atoms of Interest for Comparison:</p>"
#print  "<p><input name='ligatom' size='15' type='text' /></p>"
print " <p><input type = 'submit' value = 'Submit' /></p>"

print  " </form>"
print "</div>"

print "<div class='column'>"
print "<h2> Cofactor Substructure Based Analysis </h1>"
print "<p> Select cofactor of interest:</p>"
print  "<form action='LigPage1.py' method = 'post' target = '_blank'>"
print "<select name = 'dropdown'>"
#print "<option value='PCHILIDE'>PCHILIDE</option>"
print "        <option value='NAD'>NAD</option>"
print "        <option value='SAM'>SAM</option>"
#print "        <option value='HEM'>HEM</option>"
#print "        <option value='ATP'>ATP</option>"
print "        <option value='FAD'>FAD</option>"
print  "  </select>"

print  "  <p><input type = 'submit' value = 'Submit' /></p>"
print  " </form>"

print "</div>"


print "</div>"
print "</div>"
#print "<div class='footer'>"
#print "<p><img src='University.jpg' float='left'  width='40' height='40'<a href='https://bioenv.gu.se/english/staff?userId=xarohe'>Prof. Henrik Aronsson Group</a>, Department of Biological and Environmental Sciences, University of Gothenburg, Sweden. </p>"
#print "</div>"

print "</body>"
print "</html>"
