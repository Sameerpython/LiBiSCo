#!/usr/local/Anaconda2.7/bin/python2.7
# -*- coding: utf-8 -*-

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
    background-image: url('Title_image1_edited-111.gif') ;
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


#print "<IMG SRC='SH.jpg' ALIGN='top' width=110px height=120px >"
print "<p><b>Sameer Hassan</b></p>"
print  "<p>Postdoc Researcher</p>"
print "<p>Department of Biological & Environmental Sciences</p>"
print "<p>Email: &#032;&#115;&#097;&#109;&#101;&#101;&#114;&#046;&#104;&#097;&#115;&#115;&#097;&#110;&#064;&#107;&#105;&#046;&#115;&#101;</p>"

#print "<IMG SRC='HA.jpg' ALIGN='top' width=120px height=120px >"
print "<p><b>Mats T&#246;pel</b></p>"
print  "<p>Researcher</p>"
print "<p>Department of Marine Sciences</p>"
print "<p>Email: &#109;&#097;&#116;&#115;&#046;&#116;&#111;&#112;&#101;&#108;&#064;&#109;&#097;&#114;&#105;&#110;&#101;&#046;&#103;&#117;&#046;&#115;&#101;</p>"

#print "<IMG SRC='HA.jpg' ALIGN='top' width=120px height=120px >"
print "<p><b>Professor Henrik Aronsson</b></p>"
print  "<p>Head of Department</p>"
print "<p>Department of Biological & Environmental Sciences</p>"
print "<p>Email: &#104;&#101;&#110;&#114;&#105;&#107;&#046;&#097;&#114;&#111;&#110;&#115;&#115;&#111;&#110;&#064;&#098;&#105;&#111;&#101;&#110;&#118;&#046;&#103;&#117;&#046;&#115;&#101;</p>"

print "</div>"
print "<div class='footer'>"
print "<p><img src='University.jpg' float='left'  width='40' height='40'<a href='https://bioenv.gu.se/english/staff?userId=xarohe'>Prof. Henrik Aronsson Group</a>, Department of Biological and Environmental Sciences, University of Gothenburg, Sweden. </p>"
print "</div>"

print "</body>"
print "</html>"
