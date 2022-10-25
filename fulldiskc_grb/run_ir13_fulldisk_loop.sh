#!/bin/bash

export PATH="/home/poker/miniconda3/bin/:$PATH"

cd /home/poker/goes17/fulldiskc_grb

cp /home/ldm/data/grb-west/fulldisk/13/latest.nc /dev/shm/latest_fulldisk_13.nc
cmp /home/ldm/data/grb-west/fulldisk/13/latest.nc /dev/shm/latest_fulldisk_13.nc > /dev/null
CONDITION=$?
#echo $CONDITION

while :; do

  until [ $CONDITION -eq 1 ] ; do
#     echo same
     sleep 5
     cmp /home/ldm/data/grb-west/fulldisk/13/latest.nc /dev/shm/latest_fulldisk_13.nc > /dev/null
     CONDITION=$?
  done
#  echo different
  sleep 45
  cp /home/ldm/data/grb-west/fulldisk/13/latest.nc /dev/shm/latest_fulldisk_13.nc
  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_irc13_irc13m_and_bw.py /dev/shm/latest_fulldisk_13.nc
#  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_IR13_ircm.py /dev/shm/latest_fulldisk_13.nc
#  sleep 60
#  /home/poker/miniconda3/bin/python goes17_GRB_fulldisk_IR13_irc.py /dev/shm/latest_fulldisk_13.nc
  cmp /home/ldm/data/grb-west/fulldisk/13/latest.nc /dev/shm/latest_fulldisk_13.nc > /dev/null
  CONDITION=$?
#  echo repeat

done

