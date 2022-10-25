#!/bin/bash

export PATH="/home/poker/miniconda3/bin/:$PATH"

cd /home/poker/goes17/fulldiskc_grb

cp /home/ldm/data/grb-west/fulldisk/09/latest.nc /dev/shm/latest_fulldisk_09.nc
cmp /home/ldm/data/grb-west/fulldisk/09/latest.nc /dev/shm/latest_fulldisk_09.nc > /dev/null
CONDITION=$?
#echo $CONDITION

while :; do

  until [ $CONDITION -eq 1 ] ; do
#     echo same
     sleep 5
     cmp /home/ldm/data/grb-west/fulldisk/09/latest.nc /dev/shm/latest_fulldisk_09.nc > /dev/null
     CONDITION=$?
  done
#  echo different
  sleep 35
  cp /home/ldm/data/grb-west/fulldisk/09/latest.nc /dev/shm/latest_fulldisk_09.nc
#  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_wvc.py /dev/shm/latest_fulldisk_09.nc
  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_wvc_and_bw.py /dev/shm/latest_fulldisk_09.nc
  cmp /home/ldm/data/grb-west/fulldisk/09/latest.nc /dev/shm/latest_fulldisk_09.nc > /dev/null
  CONDITION=$?
#  echo repeat

done

