#! /bin/bash

if [ "$#" -ne 2 ]; then
    echo -e "Invalid parameters. Usage : ./sshToMachine.sh <machineType> <party>."
    exit
fi

machineType="$1"
party="$2"

lanUserName="ezpc"
lanParty0IP="40.121.40.106"
lanParty1IP="23.96.30.99"
lanParty2IP="23.96.54.155"

localUserName="ubuntu-1"
localParty0IP="10.168.236.57"
localParty1IP="10.168.236.168"
localParty2IP="10.168.236.74"

powerfulVMUserName="nishant"
powerfulVMIP="40.87.70.75"

sgxVMsUsername="mayank"
sgxVM0IP="40.118.124.169"
sgxVM1IP="13.69.21.137"
sgxVM2IP="40.68.62.177"

if [ "$machineType" == "lan" ]; then
	echo -e "Machine type given : LAN."
	if [ "$party" == "0" ];then
		echo -e "SSH to party $party. IP - $lanParty0IP."
		ssh "$lanUserName@$lanParty0IP"
	elif [ "$party" == "1" ];then
		echo -e "SSH to party $party. IP - $lanParty1IP."
		ssh "$lanUserName@$lanParty1IP"
	elif [ "$party" == "2" ];then
		echo -e "SSH to party $party. IP - $lanParty2IP."
		ssh "$lanUserName@$lanParty2IP"
	else
		echo -e "Unreognized party number - $party."
		exit
	fi
elif [ "$machineType" == "local" ]; then
	echo -e "Machine type give : Local."
	if [ "$party" == "0" ];then
		echo -e "SSH to party $party. IP - $localParty0IP."
		ssh "$localUserName@$localParty0IP"
	elif [ "$party" == "1" ];then
		echo -e "SSH to party $party. IP - $localParty1IP."
		ssh "$localUserName@$localParty1IP"
	elif [ "$party" == "2" ];then
		echo -e "SSH to party $party. IP - $localParty2IP."
		ssh "$localUserName@$localParty2IP"
	else
		echo -e "Unreognized party number - $party."
		exit
	fi
elif [ "$machineType" == "pvm" ]; then
	echo -e "Machine type give : Powerful VM."
	echo -e "SSH to IP - $powerfulVMIP."
	ssh "$powerfulVMUserName@$powerfulVMIP"
elif [ "$machineType" == "sgx" ]; then
	echo -e "Machine type give : SGX."
	if [ "$party" == "0" ];then
		echo -e "SSH to party $party. IP - $sgxVM0IP."
		ssh "$sgxVMsUsername@$sgxVM0IP"
	elif [ "$party" == "1" ];then
		echo -e "SSH to party $party. IP - $sgxV1IP."
		ssh "$sgxVMsUsername@$sgxVM1IP"
	elif [ "$party" == "2" ];then
		echo -e "SSH to party $party. IP - $sgxVM2IP."
		ssh "$sgxVMsUsername@$sgxVM2IP"
	else
		echo -e "Unreognized party number - $party."
		exit
	fi
else
	echo -e "Unrecognized machine type."
fi
