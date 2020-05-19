#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo -e "Invalid parameters. Usage : ./scpFromAllMachines.sh <machineType> <toCopyFileName> <saveFileNameSuffix>."
    exit
fi

machineType="$1"
toCopyFileName="$2"
saveFileNameSuffix="$3"

lanUserName="ezpc"
lanParty0IP="40.121.40.106"
lanParty1IP="23.96.30.99"
lanParty2IP="23.96.54.155"

localUserName="ubuntu-1"
localParty0IP="10.168.236.57"
localParty1IP="10.168.236.168"
localParty2IP="10.168.236.74"

echo -e "Machine Type : $machineType, toCopyFileName : $toCopyFileName, saveFileNameSuffix : $saveFileNameSuffix"

if [ "$machineType" == "lan" ]; then
	echo -e "Using machine type lan."
	scp "$lanUserName@$lanParty0IP:$toCopyFileName" "./tf0$saveFileNameSuffix.txt"
	scp "$lanUserName@$lanParty1IP:$toCopyFileName" "./tf1$saveFileNameSuffix.txt"
	scp "$lanUserName@$lanParty2IP:$toCopyFileName" "./tf2$saveFileNameSuffix.txt"
elif [ "$machineType" == "local" ]; then
	echo -e "Using machine type local."
	scp "$localUserName@$localParty0IP:$toCopyFileName" "./tf0$saveFileNameSuffix.txt"
	scp "$localUserName@$localParty1IP:$toCopyFileName" "./tf1$saveFileNameSuffix.txt"
	scp "$localUserName@$localParty2IP:$toCopyFileName" "./tf2$saveFileNameSuffix.txt"
else
	echo -e "Unrecognized machine type."
fi


