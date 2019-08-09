cd /Users/rachelbuttry/K2/


for i in 1 2 3 4 9 10 11 12 13 14 15 16 {21..84}
do
  echo Channel $i
  python hdf5_C8_channel_summary_1med.py $i
done

