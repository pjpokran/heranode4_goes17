#!/bin/bash

export PATH="/home/poker/miniconda3/bin/:$PATH"

cd /home/poker/goes17/fulldiskc_grb

cp /home/ldm/data/grb-west/fulldisk/14/latest.nc /dev/shm/latest_fulldisk_14.nc
cmp /home/ldm/data/grb-west/fulldisk/14/latest.nc /dev/shm/latest_fulldisk_14.nc > /dev/null
CONDITION=$?
#echo $CONDITION

while :; do

  until [ $CONDITION -eq 1 ] ; do
#     echo same
     sleep 5
     cmp /home/ldm/data/grb-west/fulldisk/14/latest.nc /dev/shm/latest_fulldisk_14.nc > /dev/null
     CONDITION=$?
  done
#  echo different
  sleep 25
  cp /home/ldm/data/grb-west/fulldisk/14/latest.nc /dev/shm/latest_fulldisk_14.nc
  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_irc_ircm_and_bw.py /dev/shm/latest_fulldisk_14.nc
#  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_IR14_ircm.py /dev/shm/latest_fulldisk_14.nc
#  sleep 60
#  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_IR14_irc.py /dev/shm/latest_fulldisk_14.nc
  cmp /home/ldm/data/grb-west/fulldisk/14/latest.nc /dev/shm/latest_fulldisk_14.nc > /dev/null
  CONDITION=$?
#  echo repeat

done

