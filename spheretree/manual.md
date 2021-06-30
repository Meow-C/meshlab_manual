# Content
* [makeTreeMedial](#makeTreeMedial)
* [makeTreeGrid](#makeTreeGrid)
* [makeTreeSpawn](#makeTreeSpawn)
* [makeTreeOctree](#makeTreeOctree)
* [makeTreeHubbard](#makeTreeHubbard)

# makeTreeMedial

./makeTreeMedial  -branch 8 -depth 1 -testerLevels 2 -numCover 10000 -minCover 5 -initSpheres 1000 -minSpheres 200 -erFact 2 -nopause -expand -merge bunny.obj

# makeTreeGrid

./makeTreeGrid  -branch 8 -depth 3 -testerLevels 2 -numCover 10000 -minCover 5 -nopause  bunny.obj

# makeTreeSpawn

./makeTreeSpawn  -branch 8 -depth 3 -testerLevels 2 -numCover 10000 -minCover 5 -nopause  bunny.obj

# makeTreeOctree

./makeTreeOctree

# makeTreeHubbard

./makeTreeHubbard  -branch 8 -depth 3 -numSamples 500 -minSamples 1 -nopause  bunny.obj
