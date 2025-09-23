#!/bin/bash

#Some cleanup
rm -rf install pcatool
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://github.com/ISSMteam/ExternalPackages/raw/refs/heads/main/pcatool.tar.gz' 'pcatool.tar.gz'

#Untar  into install
cd install 
tar -zxvf  ../pcatool.tar.gz
