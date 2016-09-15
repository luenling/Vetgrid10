import random
import argparse
import sys

def get_shuffles(replicates,num):
    num_pairs=len(replicates)/float(2)
    reps = set()
    tryouts = 0
    while ( len(reps) <= num and tryouts <= 150 ):
        tryouts +=1
        random.shuffle(replicates)
        sort_rep=sorted([ sorted(replicates[x:x+2],key=int) for x in range(0,len(replicates)-1,2)],key=lambda x: int(x[0]))
        # count proper pairs
        count_prop= 0
        for x in sort_rep:
            if (x[0] < num_pairs and x[1] >= num_pairs):
                count_prop += 1
        if count_prop > 1:
            continue
        sort_str= ",".join(["-".join(x) for x in sort_rep])
        reps.add(sort_str)
        #print sort_str
    return list(reps)

parser = argparse.ArgumentParser(description='creates a list of random shuffled comparisons for a cmh test') 

parser.add_argument("--reps1", dest="replicates1", help="the replicates to shuffle (def. \'1,2,3,4,5,6\')",default='1,2,3,4,5,6')
parser.add_argument("--reps2", dest="replicates2", help="the replicates to shuffle (def. \'7,8,9,10,11,12\')",default='7,8,9,10,11,12')
parser.add_argument("--num", dest="num", type=int, help="number of shuffles (def. 10)",default=10)

args = parser.parse_args()
replicates1 = vars(args)['replicates1'].split(',')
replicates2 = vars(args)['replicates2'].split(',')
num = vars(args)['num']
if (len(replicates1)%2 != 0 or ( len(replicates2)%2 != 0 and len(replicates2) < 1 )):
    sys.exit("not an even number of replicates")
reps1=get_shuffles(replicates1,num)
if (len(replicates2) > 1):
    reps2=get_shuffles(replicates2,num)
    reps=[ ",".join( [ reps1[x],reps2[x] ]) for x in range(0,num)]
else:
    reps=reps1
print " ".join('"' + item + '"' for item in reps)
        
    


