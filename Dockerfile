FROM ubuntu:18.04

# Update system
RUN apt-get update

# Install programs needed for the server to run
RUN apt-get install -y net-tools
RUN apt-get install -y wget 
RUN apt-get install -y sudo
RUN apt-get install -y python2.7
RUN apt-get install -y vim-tiny

RUN mkdir -p /opt/bindingdata
WORKDIR /opt/bindingdata				
RUN wget https://www.apachefriends.org/xampp-files/7.3.9/xampp-linux-x64-7.3.9-0-installer.run
RUN chmod +x xampp-linux-x64-7.3.9-0-installer.run
RUN sudo /opt/bindingdata/xampp-linux-x64-7.3.9-0-installer.run

RUN mkdir /opt/lampp/htdocs/LiBiSCo
WORKDIR /opt/lampp/htdocs/LiBiSCo
ADD httpd.conf /opt/lampp/etc/
ADD *.py /opt/lampp/htdocs/LiBiSCo/
RUN chmod +x /opt/lampp/htdocs/LiBiSCo/*.py

EXPOSE 80 
RUN /opt/lampp/xampp start
CMD sh 
# write a startup script
#RUN echo '/opt/lampp/xampp start' >> /startup.sh
#RUN echo '/usr/bin/supervisord -n' >> /startup.sh

#CMD ["sh", "/startup.sh"]
