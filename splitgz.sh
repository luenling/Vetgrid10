gunzip -c $1.gz | split -a 1 -l 200000000 - $1

