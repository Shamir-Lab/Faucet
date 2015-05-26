import tables, sets
import bloom
import strtabs # used often as string sets - keys are usually nil
import strutils, sequtils
# import nimprof
# import locks
import critbits

const
    bases = ['A', 'C', 'G', 'T']
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}.toTable()
    k = 27
    fp = 0.01
    j = 2
    read_len = 100

type
    Buff = ref object
        front, back : int
        levels : seq[StringTableRef]

proc get_kmers(r: string, k: int, kmers: var openarray[string] ) =
    # chose openarray for kmers because may want k-mers of contigs
    var i = 0
    while i < len(r)-k+1:
        kmers[i] = r[i .. i+k-1]
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
proc load_bf_sources_sinks(fname: string, numreads: int): auto =
    discard """ loads k-mers into bloom filter
        loads to sets ends of reads as potential
        sources and sinks
    """
    var
        sources = initSet[string]() # newStringTable()
        sinks = initSet[string]() # newStringTable()
        reals = initSet[string]() # newStringTable()
        bf = initialize_bloom_filter(capacity = (read_len-k+1)*numreads, error_rate = fp)   
        kmers: array[0..read_len-k, string]
        canons: array[0..read_len-k, string]
        f_hand = open(fname)
        line_no = 0
    echo(bf) 

    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo ($(line_no+1) & " read k-mers processed " & 
                $(len(sources)) & ' ' & $len(sinks) & ' ' & $len(reals))
        get_kmers(line,k,kmers)
        get_canons(kmers,canons)
        sources.incl(canons[0]) # ] = nil
        sinks.incl(canons[read_len-k]) # ] = nil
        for i,value in @canons:
            if value!=nil:
                bf.insert(value)
                # reals only for debugging:
                reals.incl(value) # ]=nil
        inc(line_no)
    f_hand.close()
    echo($len(reals) & " real k-mers loaded")
    (bf,sources,sinks,reals)


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
    front.incl(kmers[0])
    for i in 0..j:
        advance_front(i, kmers, front, bf)


proc get_candidate_paths(filename: string, bf: object; rc=false): auto =

    discard """ scan reads to find candidate j+1 length paths. 
        finds nodes having descendents at level j+1 that differ from read sequence
        used to check against reals later, to know which are false joins, vs. true
        alternate paths; returns list of candidate lists corr. to all candidate per read
    """
    var
        f_hand = open(filename)
        cands : CritBitTree[void] # newStringTable()
        line_no = 0
        read: string
        kmers: array[0..read_len-k, string]
        front : CritBitTree[void]
        front_prefs : CritBitTree[void]
        next_real = ""
        front_real = ""
        back_suffix = ""

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

        for ind, kmer in @kmers:
            # can only set when next real is known
            # otherwise every front node is alt path
            if ind < read_len-k:
                next_real = kmers[ind+1]
                back_suffix = next_real[j..k-1]
            # echo(front)
            for s in front.items:
                # below only holds for alts
                # echo("back suffix: " & back_suffix & " front prefix: " & s[0..k-j-1])
                if s[0..k-j-1]!=back_suffix:
                    # alt path is source k-mer concat with last j+1 chars
                    # of node in front
                    # echo("position " & $ind)
                    cands.incl(kmer) # & s[k-j-1..k-1])
                    front.excl(s)
            # echo($(ind+j+1))
            advance_front(ind+j+1, kmers, front, bf)

            # test if alt paths at front
            # when found, are added to cands, removed from front
            # advance front

        # echo($len(cands) & " cands: "& $cands)
        inc(line_no)
        
    if rc==true:
        echo("got rev. com. candidates")
    else:
        echo("got forward candidates") 
    f_hand.close()
    return cands


when isMainModule:
    var 
        reads_file = "/vol/scratch/rozovr/chr20.c10.reads.100k"
        #reads_file = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.100k"
        (bf,sources,sinks,reals)=load_bf_sources_sinks(reads_file, 100_000)
        bf_cands = get_candidate_paths(reads_file, bf)

    # solution with threads
    # load BF using 1 thread
    # find candidates with 10 threads - 
    # each thread reads from block of reads (based on index parameter)
    # all can query (read only from) bf
    # each returns candidate set, take union at end
