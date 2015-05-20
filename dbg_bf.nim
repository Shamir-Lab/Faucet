import tables, sets
import bloom
import strtabs # used often as string sets - keys are usually nil
import strutils, sequtils
import nimprof

const
    bases = ['A', 'C', 'G', 'T']
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}.toTable()
    k = 27
    fp = 0.01
    j = 10
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
        sources = newStringTable()
        sinks = newStringTable()
        reals = newStringTable()
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
        sources[canons[0]] = nil
        sinks[canons[read_len-k]] = nil
        for i,value in @canons:
            if value!=nil:
                bf.insert(value)
                # reals only for debugging:
                reals[value]=nil
        inc(line_no)
    f_hand.close()
    echo($len(reals) & " reals k-mers loaded")
    (bf,sources,sinks,reals)


proc get_empty_buff(j: int): Buff = 
    # newSeqWith seen at http://forum.nim-lang.org/t/1161
    var levels : seq[StringTableRef] = newSeqWith(j+2, newStringTable())
    # init with one extra empty level, that's rotated from front to back
    # to avoid allocations
    return Buff(front:j,back:0,levels:levels)

proc print_buff_info(buff: Buff) =
    for i in 0..len(buff.levels)-1:
        var x = 0
        echo("buff level " & $i & ": ")
        for key in buff.levels[i].keys:
            echo($(x+1) & ' ' & key)
            inc(x)

proc init_read_buff(source: string, buff: var Buff, bf: object) = 
    var 
        roots, next : StringTableRef
        test_kmer, canon : string
    buff.front = j
    buff.back = 0
    # buff.levels[0] = newStringTable()
    buff.levels[0][source]=nil
    # buff.levels[1] = newStringTable()
    for level in 0..j:
        roots = buff.levels[level]
        next = buff.levels[level+1]
        if next == nil:
            break
        for root in roots.keys:
            for b in bases: 
                test_kmer = root[1..k-1] & b
                canon = min(test_kmer, get_rc(test_kmer))
                if bf.lookup(canon)==true:
                    next[test_kmer]=nil
        buff.levels[level+1] = next
        # print_buff_info(buff)

proc get_buffer_level(buff: Buff, j,level: int): TableRef[string,seq[string]] =
    discard """ gets buffer contents from chosen level (in [0,j])
        returns dictionary containing invariant k-j as keys
        and k-mers including them as values - e.g., level = 0
        contains (k-j)-mers as suffixes, level = j+1 contains 
        them as prefixes  
    """
    var inv: string
    result = newTable[string,seq[string]]()

    for kmer in buff.levels[level].keys:
        if level != 0:
            inv = kmer[j - level .. ^(level+1)]
            # inv = kmer[j+level .. k-level-1]
        else:
            inv = kmer[j .. k-1]

        # echo("inv in get_level " & inv)
        if hasKey(result, inv):
            result.mget(inv).add(kmer)
        else:
            result[inv] = @[kmer]

proc get_alt_paths_from_buff(comms: StringTableRef, next_real: string,
 backs, fronts: TableRef[string,seq[string]], buff: Buff): auto =
    discard """ given (k-j)-mers of alt paths, gets their start and end
        k-mers to create paths to check (list returned)
    """
    var
        pref, path : string
        ends = newSeq[string]()
    result = newSeq[string]() # newStringTable()
    # echo($len(comms) & ' ' & next_real)
    for comm in comms.keys:
        if comm != next_real:
            # echo("\n" & comm)
            # echo("backs")
            # echo(backs[comm])
            # echo("\n" & $backs)
            # echo("fronts")
            # echo(fronts[comm])
            # echo("\n" & $fronts)
            pref = backs[comm][0]
            ends = fronts[comm]
            for fr in ends:
                path = pref & fr[k-(1+j) .. k-1]
                result.add(path) # [path]=nil


proc clean_front(buff: var Buff, fronts:TableRef[string,seq[string]], 
    comms:StringTableRef) =
    discard """ removes all nodes starting with alt. (k-j)-mers
        from buffer front
    """
    var new_front = newStringTable()
    for fr in fronts.keys:
        # instead of discarding, only insert non-alts
        # then replace buffer front
        if (not hasKey(comms, fr)):
            for key in fronts[fr]:
                new_front[key]=nil
    buff.levels[buff.front]=new_front

proc advance_buffer(buff: var Buff, bf: object) = 
    discard """ does BFS using Bloom filter bf from loaded buffer
        extends each front node one step, removes first level
    """
    var 
        fronts = newStringTable()
        next = newStringTable()
        test_kmer, canon: string
    fronts = buff.levels[j]
    for node in fronts.keys:
        for b in bases:
            test_kmer = node[1 .. k-1] & b
            # note: python code was missing next line
            canon = min(test_kmer, get_rc(test_kmer))
            if bf.lookup(canon)==true:
                next[test_kmer]=nil
    buff.front = (buff.front+1) mod (j+2)
    # echo("new front is " & $buff.front)
    buff.levels[buff.front] = next
    buff.back = (buff.back+1) mod (j+2)
    # echo("new back is " & $buff.back)


proc get_candidate_paths(filename: string, bf: object; rc=false): auto =

    discard """ scan reads to find candidate j+1 length paths. 
        finds nodes having descendents at level j+1 that differ from read sequence
        used to check against reals later, to know which are false joins, vs. true
        alternate paths; returns list of candidate lists corr. to all candidate per read
    """
    var
        f_hand = open(filename)
        cands = newStringTable()
        line_no = 0
        read: string
        kmers: array[0..read_len-k, string] # remember both ends included in ranges
        buff = get_empty_buff(j)
        backs : ref Table[string,seq[string]]
        fronts : ref Table[string,seq[string]]
        comms = newStringTable()
        # alts : newStringTable()
        next_real: string
        alt_paths = newSeq[string]()

    echo("in cand paths")
    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo($(line_no+1) & ' ' & $len(cands))
        read = strip(line)
        if rc:
            read = get_rc(read)
        get_kmers(read, k, kmers)
        # echo(read)
        init_read_buff(kmers[0], buff, bf)
        # print_buff_info(buff)
        for ind, value in @kmers:
            backs = get_buffer_level(buff,j,buff.back)
            # echo("backs " & $len(backs))
            # for key in backs.keys:
            #         echo(key & ": " & backs[key])
            
            fronts = get_buffer_level(buff,j,buff.front)
            # echo("fronts" & $len(fronts))
            # for key in fronts.keys:
            #     echo(key & ": " & fronts[key])
            
            for key in backs.keys: # get intersection
                if haskey(fronts, key):
                    # echo(key & ": " & backs[key])
                    # echo(key & ": " & fronts[key])
                    comms[key]=nil
                    # echo(comms)
            if len(comms) >= 2:
                if ind < read_len-k-1:
                    next_real = kmers[ind+1][j..k-1]
                else:
                    next_real = ""
                # note next_real not removed from string table
                # comms because I don't know how...
           
                alt_paths = get_alt_paths_from_buff(comms, next_real, backs, fronts, buff)
                if len(alt_paths) > 0:
                    # echo(alt_paths)
                    for alt in alt_paths:
                        cands[value[0] & alt] = nil
                    clean_front(buff,fronts,comms)
            advance_buffer(buff,bf)
            comms = newStringTable()

        inc(line_no)
        # clear out buffer state before re-use
        for i in 0 .. len(buff.levels)-1:
            buff.levels[i] = newStringTable() 
            #### thought would be faster, but doesn't work:
            #### clear(buff.levels[i], modeCaseSensitive)
    if rc==true:
        echo("got rev. com. candidates")
    else:
        echo("got forward candidates") 
    f_hand.close()
    return cands

# def get_candidate_paths(filename,bf,rc=False):
#     """ scan reads to find candidate j+1 length paths. 
#         finds nodes having descendents at level j+1 that differ from read sequence
#         used to check against reals later, to know which are false joins, vs. true
#         alternate paths; returns list of candidate lists corr. to all candidate per read
#     """
#     cands = []
#     line_no = 0
#     with open(filename) as f:
#         for line in f:
#             if (line_no+1)%10000==0:
#                 print line_no+1, len(cands)
#             read = line.rstrip()
#             if rc:
#                 kmers = get_kmers(get_rc(read),k)
#             else:
#                 kmers = get_kmers(read,k)
#             buff = get_j_forward_buff(kmers[0],bf,j)
            
#             for ind, kmer in enumerate(kmers):
#                 backs = get_buffer_level(buff,j,0) 
#                 fronts = get_buffer_level(buff,j,j)
#                 comms = (set(backs.keys())).intersection(set(fronts.keys()))

#                 if len(comms) >= 2: 
#                     if ind < read_len-k:
#                         next_real = kmers[ind+1][j:] # k-j suffix from next real k-mer  
#                     else: # if read's end reached, don't know what next real is
#                         next_real = set()                   
#                     alts = comms - set([next_real]) 
#                     alt_paths = get_alt_paths_from_buff(alts, backs, fronts, buff)
#                     if alt_paths:
#                         for alt in alt_paths:
#                             alt = kmer[0] + alt # add first letter of previous k-mer
#                         cands.append(alt_paths)
#                         clean_front(buff,fronts,alts)

#                 advance_buffer(buff,bf)
#                 if ind < read_len-k: # need to think more about read ends
#                     clean_back(buff,kmers[ind+1])   
#             line_no +=1
#     if rc:
#         print "got rc candidates"
#     else:
#         print "got forward candidates"
#     return cands


when isMainModule:
    var 
        reads_file = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.100k"
        (bf,sources,sinks,reals)=load_bf_sources_sinks(reads_file, 100_000)
        bf_cands = get_candidate_paths(reads_file, bf)
    # var read = "ACGTTCGTTTGACACTTCGTTTGTCGTTTGGTTCGTTGTTCGTT"
    # echo reverse(read)
    # echo get_rc(read)
    # echo read # original isn't changed
        
