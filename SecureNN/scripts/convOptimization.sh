#!/bin/bash

MinFilterSize=1
MaxFilterSize=15
BatchSize=64
PartyNum=$1
cd ..
echo -e "PartyNum = $PartyNum"
for ((curFilterSize=$MinFilterSize;curFilterSize<=$MaxFilterSize;curFilterSize++)); do
	echo -e "Starting new FilterSize. curFilterSize = $curFilterSize."
	if [ "$PartyNum" == "0" ]; then
		./SecureNN.out 3PC 0 files/parties_LAN files/keyA files/keyAB "$BatchSize" "$curFilterSize"
	elif [ "$PartyNum" == "1" ]; then
		./SecureNN.out 3PC 1 files/parties_LAN files/keyB files/keyAB "$BatchSize" "$curFilterSize"
	else
		./SecureNN.out 3PC 2 files/parties_LAN files/keyB files/keyAB "$BatchSize" "$curFilterSize"
	fi
done
