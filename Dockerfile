FROM ubuntu:18.04

# Update system
RUN apt-get update

# Install programs needed for the server to run
RUN apt-get install -y net-tools
RUN apt-get install -y wget 
RUN apt-get install -y sudo
RUN apt-get install -y vim-tiny
RUN apt-get install -y libgl1-mesa-glx

RUN mkdir -p /opt/bindingdata
WORKDIR /opt/bindingdata

# Install XAMPP
RUN wget https://www.apachefriends.org/xampp-files/7.3.9/xampp-linux-x64-7.3.9-0-installer.run
RUN chmod +x xampp-linux-x64-7.3.9-0-installer.run
RUN sudo /opt/bindingdata/xampp-linux-x64-7.3.9-0-installer.run

# Install Anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh
RUN chmod +x Anaconda2-2019.10-Linux-x86_64.sh
RUN sudo ./Anaconda2-2019.10-Linux-x86_64.sh -b -p /usr/local/Anaconda2.7
RUN echo "PATH=/usr/local/Anaconda2.7/bin:$PATH" >> /etc/profile

# Install python dependencies
RUN /usr/local/Anaconda2.7/bin/conda install -y -c conda-forge biopython

RUN export PATH=/usr/local/Anaconda2.7/bin:$PATH

RUN mkdir /opt/lampp/htdocs/LiBiSCo
WORKDIR /opt/lampp/htdocs/LiBiSCo
ADD httpd.conf /opt/lampp/etc/
ADD *.py /opt/lampp/htdocs/LiBiSCo/
ADD *.gif /opt/lampp/htdocs/LiBiSCo/
ADD *.png /opt/lampp/htdocs/LiBiSCo/
RUN chmod +x /opt/lampp/htdocs/LiBiSCo/*.py

# RUN echo -e "[Unit]\nDescription=XAMPP\n\n[Service]\nExecStart=/opt/lampp/lampp start\nExecStop=/opt/lampp/lampp stop\nType=forking\n\n[Install]\nWantedBy=multi-user.target"

RUN echo -e "[Unit]\nDescription=XAMPP\n\n[Service] \n\
ExecStart=/opt/lampp/lampp start\nExecStop=/opt/lampp/lampp stop \n\
Type=forking\n\n[Install]\nWantedBy=multi-user.target" > /etc/systemd/system/xampp.service 

EXPOSE 80 
CMD /opt/lampp/xampp start; sh 
# write a startup script
#RUN echo '/opt/lampp/xampp start' >> /startup.sh
#RUN echo '/usr/bin/supervisord -n' >> /startup.sh

#CMD ["sh", "/startup.sh"]
