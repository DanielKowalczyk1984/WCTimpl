
for n in 100 ; do
for m in  3   ; do
sed -i.bak 's/;/,/g' WCT_AHV_$m\_$n.csv
done
done
