$SEED=73073

./build/iqtree3 -s test_scripts/test_data/small.fa -p test_scripts/test_data/small.nex -m "MFP+MERGE" -T 1 -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -B 1000 -T 1 -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.nex -B 1000 -T 1 -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.nex -B 1000 -T 1 -m MFP+MERGE -rcluster 10 --prefix test_scripts/test_data/turtle.merge -seed $SEED

Get-Content test_scripts/test_data/turtle.fa.treefile, test_scripts/test_data/turtle.nex.treefile |
    Set-Content test_scripts/test_data/turtle.trees
./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.merge.best_scheme.nex -z test_scripts/test_data/turtle.trees -zb 10000 -au -n 0 --prefix test_scripts/test_data/turtle.test -seed $SEED -T 1

./build/iqtree3 -s test_scripts/test_data/turtle.fa -m GTR+F+I+R3+T -te test_scripts/test_data/turtle.trees -T 1 --prefix test_scripts/test_data/turtle.mix -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.nex.best_scheme.nex -z test_scripts/test_data/turtle.trees -n 0 -wpl --prefix test_scripts/test_data/turtle.wpl -seed $SEED -T 1

./build/iqtree3 -s test_scripts/test_data/turtle.fa -S test_scripts/test_data/turtle.nex --prefix test_scripts/test_data/turtle.loci -T 1 -seed $SEED

./build/iqtree3 -t test_scripts/test_data/turtle.nex.treefile --gcf test_scripts/test_data/turtle.loci.treefile -s test_scripts/test_data/turtle.fa --scf 100 -seed $SEED -T 1

./build/iqtree3 -t test_scripts/test_data/turtle.fa.treefile --gcf test_scripts/test_data/turtle.loci.treefile -s test_scripts/test_data/turtle.fa --scf 100 -seed $SEED -T 1

# link-exchange-rates model

./build/iqtree3 -s test_scripts/test_data/turtle.fa -m "MIX{GTR+FO,GTR+FO}" --link-exchange-rates --prefix test_scripts/test_data/turtle.mix.link -seed $SEED -T 1

./build/iqtree3 -s test_scripts/test_data/turtle.fa -m "MIX{GTR{1,1,1,1,1,1}+FO,GTR{1,1,1,1,1,1}+FO}" --link-exchange-rates --prefix test_scripts/test_data/turtle.mix.jc.link -seed $SEED -T 1

./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.nex -g test_scripts/test_data/turtle.constr.tree --prefix test_scripts/test_data/turtle.nex.constr -T 1 -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -p test_scripts/test_data/turtle.nex -g test_scripts/test_data/turtle.constr.tree2 -B 1000 -alrt 1000 --prefix test_scripts/test_data/turtle.nex.constr2 -T 1 -seed $SEED

./build/iqtree3 -s test_scripts/test_data/turtle.fa -m "MIX+MF" --prefix test_scripts/test_data/turtle.mixfinder -T 1 -seed $SEED
