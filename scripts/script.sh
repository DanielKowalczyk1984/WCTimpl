
for m in   3; do
for n in  100 ; do
for file in ../$n*$m/*.txt;do
./wct  -p 1 -l 900 -S 0 -z 1 -b 1 -Z 1 -f 10 -c 0 $file $m;
echo $m;
done
done
done
