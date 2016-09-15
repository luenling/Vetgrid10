import random
import itertools
case=["1","2","3"]
case_tests=[]
for i in itertools.permutations(case,2):
    case_tests.append(":".join(i))
base=["4","5","6"]
base_tests=[]
for i in itertools.permutations(base,2):
    base_tests.append(":".join(i))

all = base_tests + case_tests
n=0
count = 0
all_set=set([])
while (n <= 10) and (count < 100):
    all_set.add(",".join(sorted(random.sample(all,3))))
    n = len(all_set)
    count += 1

