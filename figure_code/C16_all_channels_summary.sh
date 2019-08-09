cd /Users/rachelbuttry/K2/


for i in 1 2 3 4 13 14 15 16 {21..84}
do
  echo Channel $i
  python hdf5_C16_channel_summary_1med.py $i
done

