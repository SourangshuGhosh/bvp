#! /bin/bash
#
python3 bvp_01.py > bvp_01.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
python3 bvp_02.py > bvp_02.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
python3 bvp_03.py > bvp_03.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
python3 bvp_04.py > bvp_04.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
python3 bvp_05.py > bvp_05.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
