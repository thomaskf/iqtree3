$SEED=73073

# ./build/iqtree3 -s test_scripts/test_data/small.fa -p test_scripts/test_data/small.nex -m "MFP+MERGE" -T 1 -seed $SEED

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

## amino acid test cases
# the data set is a subset of the turtle data set
echo "Running amino acid test cases..."
AA_FASTA=test_scripts/test_data/turtle_aa.fasta
AA_NEX=test_scripts/test_data/turtle_aa.nex
AA_prefix=test_scripts/test_data/turtle_aa

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -B 1000 -T 1 -seed $SEED

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p $AA_NEX -B 1000 -T 1 -seed $SEED

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p $AA_NEX -B 1000 -T 1 -m MFP+MERGE -rcluster 10 --prefix test_scripts/test_data/turtle_aa.merge -seed $SEED

Get-Content $AA_FASTA.treefile, $AA_NEX.treefile |
    Set-Content test_scripts/test_data/turtle_aa.trees
"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p test_scripts/test_data/turtle_aa.merge.best_scheme.nex -z test_scripts/test_data/turtle_aa.trees -zb 10000 -au -n 0 --prefix test_scripts/test_data/turtle_aa.test -seed $SEED -T 1

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -m GTR+F+I+R3+T -te test_scripts/test_data/turtle_aa.trees -T 1 --prefix test_scripts/test_data/turtle_aa.mix -seed $SEED

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p $AA_NEX.best_scheme.nex -z test_scripts/test_data/turtle_aa.trees -n 0 -wpl --prefix test_scripts/test_data/turtle_aa.wpl -seed $SEED -T 1

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -S $AA_NEX --prefix test_scripts/test_data/turtle_aa.loci -T 1 -seed $SEED

"$BUILD_DIR/iqtree3"  -t $AA_NEX.treefile --gcf test_scripts/test_data/turtle_aa.loci.treefile -s $AA_FASTA --scf 100 -seed $SEED -T 1

"$BUILD_DIR/iqtree3"  -t $AA_FASTA.treefile --gcf test_scripts/test_data/turtle_aa.loci.treefile -s $AA_FASTA --scf 100 -seed $SEED -T 1

# link-exchange-rates model

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -m "MIX{LG+F,LG+F}" --link-exchange-rates --prefix test_scripts/test_data/turtle_aa.mix.link -seed $SEED -T 1

#"$BUILD_DIR/iqtree3"  -s $AA_FASTA -m "MIX{GTR{1,1,1,1,1,1}+FO,GTR{1,1,1,1,1,1}+FO}" --link-exchange-rates --prefix test_scripts/test_data/turtle_aa.mix.jc.link -seed $SEED -T 1

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p $AA_NEX -g test_scripts/test_data/turtle.constr.tree --prefix $AA_NEX.constr -T 1 -seed $SEED

"$BUILD_DIR/iqtree3"  -s $AA_FASTA -p $AA_NEX -g test_scripts/test_data/turtle.constr.tree2 -B 1000 -alrt 1000 --prefix $AA_NEX.constr2 -T 1 -seed $SEED

#"$BUILD_DIR/iqtree3"  -s $AA_FASTA -m "MIX+MF" --prefix test_scripts/test_data/turtle_aa.mixfinder -T 1 -seed $SEED # mixture finder is not supported for amino acid data
