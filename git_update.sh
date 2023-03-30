#!/bin/sh

#init git 
#git init

#add files 
git add *.R *.sh 

#commit files 
git commit -m "update git update script!"

#set the origin, only for the first time
#git remote add origin git@github.com:ychenatsinca/WuF_site.git
#add branch name, here is main 
#git branch -M main
#push commit files to the server/origin as master or branch 
git push -u origin main


