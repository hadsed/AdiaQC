#!/bin/bash
cd ../../

for a in {1..1000}
do
    echo -e "\n\nStarting iteration $a.\n\n"
    python2 run.py -p hopfield/tetris_hebbrule
    python2 run.py -p hopfield/tetris_storkeyrule
    python2 run.py -p hopfield/tetris_projectionrule
done
