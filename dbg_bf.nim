import tables, sets
import bloom
from strutils import strip
# import nimprof
# import locks
import critbits

const
    bases = ['A', 'C', 'G', 'T']
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}.toTable()
    k = 27
    fp = 0.01
    j = 1
    read_len = 100

proc get_kmers(r: string, sub_len: int, kmers: var openarray[string] ) =
    # chose openarray for kmers because may want k-mers of contigs
    var i = 0
    while i < len(r)-sub_len+1:
        kmers[i] = r[i .. i+sub_len-1]
        i = i+1

proc reverse(s: string): string =
    # from http://goran.krampe.se/2014/10/29/nim-and-oo/
    result = s
    var x = 0
    var y = s.high
    while x < y:
        swap(result[x], result[y])
        dec(y)
        inc(x)

proc get_rc(dna: string): string = 
    result = reverse(dna) 
    for i in 0 .. dna.high:
        result[i] = complements[result[i]]

proc get_canons(kmers: openarray[string], canons: var openarray[string]) =
    for i,value in @kmers:
        canons[i]=min(value,get_rc(value))


#### note use of auto to determine types inside tuple returned
proc load_bf(fname: string, sub_len, numreads: int): auto =
    discard """ loads k-mers into bloom filter
    """
    var
        bf = initialize_bloom_filter(capacity = (read_len-sub_len+1)*numreads, error_rate = fp)   
        kmers: seq[string] = newSeq[string](read_len-sub_len+1)
        canons: seq[string] = newSeq[string](read_len-sub_len+1)

        f_hand = open(fname)
        line_no = 0
    
    echo(bf) 

    for line in f_hand.lines:
        if (line_no + 1) mod 30_000==0:
            echo ($(line_no+1) & " read k-mers inserted to BF ")
        get_kmers(line,sub_len,kmers)
        get_canons(kmers,canons)
        for i,value in @canons:
            if value!=nil:
                bf.insert(value)
                
        inc(line_no)
    f_hand.close()
    return bf 
    
proc load_reals(fname: string, sub_len: int): auto =
    discard """ loads k-mers into set
    """
    var
        reals: CritBitTree[void]
        kmers: seq[string] = newSeq[string](read_len-sub_len+1)
        canons: seq[string] = newSeq[string](read_len-sub_len+1)

        f_hand = open(fname)
        line_no = 0
    
    for line in f_hand.lines:
        if (line_no + 1) mod 30_000==0:    
            echo ($(line_no+1) & " read k-mers inserted to reals set")

        get_kmers(line,sub_len,kmers)
        get_canons(kmers,canons)
        for i,value in @canons:
            if value!=nil:    
                reals.incl(value)
        inc(line_no)
    f_hand.close()
    return reals 


proc load_alt_extensions(curr_real, next_real: string, 
    bf: object, next_set: var CritBitTree[void]) = 
    var 
        test_kmer = curr_real[1..k-1]
        canon : string
    # echo("curr_real: " & curr_real & " next_real: " & next_real)
    if test_kmer == "": return

    for b in bases:
        if b == 'A': 
            test_kmer.add(b)
        else:
            test_kmer[k-1]=b     

        # only want to consider alternatives to read sequence
        if test_kmer == next_real: continue

        canon = min(test_kmer, get_rc(test_kmer))
        if bf.lookup(canon)==true:
            next_set.incl(test_kmer)

proc advance_front(read_pos: int, kmers: openarray[string],
 front: var CritBitTree[void], bf: object) = 
    var 
        next : CritBitTree[void]
        front_kmer = ""
        next_real = ""
    if read_pos < read_len - k:
        front_kmer = kmers[read_pos]
    if read_pos < read_len - k - 1:
        next_real = kmers[read_pos+1]
    # echo(front_kmer & " " & next_real)
    for s in front.items:
        # alt extensions --> diff from next (read) kmer
        # put alt extensions of curr front nodes 
        # diff from next_real into next
        load_alt_extensions(s,next_real,bf,next)

    # put alt extensions of real at front into next
    load_alt_extensions(front_kmer,next_real,bf,next)
    front = next
    # echo(front)

proc load_front(kmers: openarray[string], front: var CritBitTree[void], bf: object) = 
    if len(front)>0:
        for s in front.items:
            front.excl(s)
    front.incl(kmers[0])
    for i in 0..j:
        advance_front(i, kmers, front, bf)


proc get_candidate_paths(filename: string, bf: object; rc=false): auto =

    discard """ scan reads to find candidate j+1 length paths. 
        finds nodes having descendents at level j+1 that differ from read sequence
        used to check against reals later, to know which are false joins, vs. true
        alternate paths; returns set of candidate (k+j+1) length paths
    """
    var
        f_hand = open(filename)
        cands : CritBitTree[void] # newStringTable()
        line_no = 0
        # cand_cnt = 0
        read: string
        kmers: array[0..read_len-k, string]
        front : CritBitTree[void]
        next_real = ""
        back_suffix = ""
        first_buff_pref = ""
        added : bool

    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo($(line_no+1) & ' ' & $len(cands))
        read = strip(line)
        if rc:
            read = get_rc(read)
        get_kmers(read, k, kmers)
        # echo(read)
        load_front(kmers, front, bf)
        # if line_no+1==10: break
        added = false

        for ind, kmer in @kmers:
            # can only set when next real is known
            # otherwise every front node is alt path
            if ind < read_len-k:
                next_real = kmers[ind+1]
                back_suffix = next_real[j..k-1]
            first_buff_pref = kmer[1..k-1]
            # echo(front)
            for s in front.items:
                # below only holds for alts
                # echo("back suffix: " & back_suffix & " front prefix: " & s[0..k-j-1])
                if s[0..k-j-1]!=back_suffix:
                    if not added:
                        cands.incl(next_real)
                        added = true
                    front.excl(s)

            #         # alt path is source k-mer concat with last j+1 chars
            #         # of node in front
            #         # echo("position " & $ind)
            #         cands.incl(first_buff_pref & s[k-j-1]) # candidate to check is first k-mer in buffer
            #         inc(cand_cnt)
            #         front.excl(s)
            # # echo($(ind+j+1))
            advance_front(ind+j+1, kmers, front, bf)

        # echo($len(cands) & " cands: "& $cands)
        inc(line_no)
        
    if rc==true:
        echo("got rev. com. candidates")
    else:
        echo("got forward candidates") 
    f_hand.close()
    return cands #, cand_cnt)


when isMainModule:
    var 
        reads_file = "/vol/scratch/rozovr/chr20.c30.orc_out.reads.tail.1M" #"/vol/scratch/rozovr/chr20.c10.reads.1M"
        #reads_file = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.100k"
        bf1 = load_bf(reads_file, k, 1_000_000) # ,sources,sinks,reals)
        # reals = load_reals(reads_file, k)
        # bf2 = load_bf(reads_file, k+j+1, 1_000_000) # ,sources,sinks,reals)

        bf_cands = get_candidate_paths(reads_file, bf1)
        # num_fps = 0
    # echo("cands paths list size: " & $cnd_cnt & ", cands juncs set size: " & $len(bf_cands))
    echo("cands real junc set size: " & $len(bf_cands))
    # for s in bf_cands.items:
    #     if not reals.contains(s) and not reals.contains(get_rc(s)):
    #         inc(num_fps)

    # echo($num_fps & " false junctions found") 

        # rc_cands = get_candidate_paths(reads_file, bf1, true)

    # solution with threads
    # load BF using 1 thread
    # find candidates with 10 threads - 
    # each thread reads from block of reads (based on index parameter)
    # all can query (read only from) bf
    # each returns candidate set, take union at end
