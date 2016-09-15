#!/bin/bash
PATH=/Volumes/Temp/Lukas/Tools/bin:/Library/Frameworks/Python.framework/Versions/2.7/bin:/Library/Frameworks/Python.framework/Versions/2.7/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin
FILENAME="bam_file_database_"`date +"_%d.%m.%y"`
DLDIR=/Users/vetgrid10/Public/BamStorage
cd $DLDIR
google docs get -n "BAM file database" --format=xls --dest=$FILENAME
exit

