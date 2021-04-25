#!/usr/bin/python

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
    background-image: url('Pic.gif') ;
    background-size: cover;
    padding-left: 80px;
    
}
}
"""
print """
img {
    float: left;
    
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
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href='Contact.py'>Contact</a></li>"
print "</ul>"
#Body size Determination


print "<div align='center'>"
print "<h1> Please send any comments or suggestions about the LiBiSCo server to:</h1>"
#print "</br>"
print "<div id='container1'>"

#print "<h2> Prof. Henrik Aronsson - henrik.aronsson@bioenv.gu.se </h2>"
#print "<p><img  src='HA.jpg'  style='width:120px;height:120px;margin-right:15px;'>"
#print "<b>Prof. Henrik Aronsson</b> - henrik.aronsson@bioenv.gu.se</p>" 

print "<IMG SRC='HA.jpg' ALIGN='top' width=120px height=120px >"
print "<p><b>Professor Henrik Aronsson</b></p>"
print  "<p>Head of Department</p>"
print "<p>Department of Biological & Environmental Sciences</p>"
print "<p>Email: henrik.aronsson@bioenv.gu.se</p>"



print "<div class='row'>"
print "  <div class='column'>"
print "<IMG SRC='SH.jpg' ALIGN='top' width=110px height=120px >"
print "<p><b>Sameer Hassan</b></p>"
print  "<p>Postdoc Researcher</p>"
print "<p>Department of Biosciences and Nutrition</p>"
print "<p>Email: sameer.hassan@ki.se</p>"
print"</div>"
print"<div class='column'>"
print "<IMG SRC='MT.png' ALIGN='top' width=110px height=120px >"
print "<p><b>Mats Topel</b></p>"
print  "<p>Researcher</p>"
print "<p>Department of Marine Sciences</p>"
print "<p>Email: mats.topel@marine.gu.se</p>"
print "</div>"
print "</div>"
#print "Professor , Head of Department ", "<br/"
#print "<h2>Sameer Hassan - sameer.hassan@bioenv.gu.se </h2>"
print "</div>"

print "</body>"
print "</html>"
