import tables
import bloom
import strtabs # used often as string sets - keys are usually nil
import strutils, sequtils

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
        kmers: array[0..read_len-k+1, string]
        canons: array[0..read_len-k+1, string]
        f_hand = open(fname)
        line_no = 0
    echo(bf) 

    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo ($(line_no+1) & " read k-mers processed " & 
                $(len(sources)) & " " & $len(sinks) & " " & $len(reals))
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
    (bf,sources,sinks,reals)


proc get_empty_buff(j: int): Buff = 
    # newSeqWith seen at http://forum.nim-lang.org/t/1161
    var levels : seq[StringTableRef] = newSeqWith(j+1, newStringTable())
    return Buff(front:0,back:0,levels:levels)

proc init_read_buff(source: string, buff: var Buff, bf: object) = 
    var 
        roots, next : StringTableRef
        test_kmer, canon : string

    buff.levels[0][source]=nil
    for level in 0..j:
        roots = buff.levels[level]
        next = buff.levels[level+1]
        if next == nil:
            break
        for root in roots.keys:
            for b in bases: 
                test_kmer = root[1..root.len] & b
                canon = min(test_kmer, get_rc(test_kmer))
                if bf.lookup(canon)==true:
                    next[test_kmer]=nil

proc print_buff_info(buff: Buff) =
    for i in 0..len(buff.levels)-1:
        if buff.levels[i] == nil:
            echo "empty buffer level"
        else:
            # echo($i & ": " & $len(buff.levels[i]))
            var x = 0
            echo($i & ": ")
            for key in buff.levels[i].keys:
                echo($(x+1) & " " & key)
                inc(x)

proc get_buffer_level(buff: Buff, j,level: int): TableRef[string,seq[string]] =
    discard """ gets buffer contents from chosen level (in [0,j])
        returns dictionary containing invariant k-j as keys
        and k-mers including them as values - e.g., level = 0
        contains (k-j)-mers as suffixes, level = j+1 contains 
        them as prefixes  
    """
    var inv: string
    result[] = initTable[string,seq[string]]()
    for kmer in buff.levels[level].keys:
        # if level != 0:
        inv = kmer[j+level .. k-level-1]

        if hasKey(result, inv):
            result.mget(inv).add(kmer)
        else:
            result[inv] = @[kmer]
    
# def get_buffer_level(buff, j, level):
    # """ gets buffer contents from chosen level (in [0,j])
    #     returns dictionary containing invariant k-j as keys
    #     and k-mers including them as values - e.g., level = 0
    #     contains (k-j)-mers as suffixes, level = j+1 contains 
    #     them as prefixes  
    # """
#     kmers = list(buff[level])
#     invars = {}
#     for kmer in kmers:
#         if level != 0:
#             inv = kmer[j - level : -level]
#         else:
#             inv = kmer[j - level :]

#         if inv in invars:
#             invars[inv].append(kmer)
#         else:
#             invars[inv]=[kmer]
#     return invars

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
        kmers: array[0..read_len-k+1, string]
        buff = get_empty_buff(j)
        backs = newTable[string,seq[string]]()
        fronts = newTable[string,seq[string]]()

    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo($(line_no+1) & " " & $len(cands))
        read = strip(line)
        if rc:
            read = get_rc(read)
        get_kmers(read, k, kmers)
        init_read_buff(kmers[0], buff, bf)
        # print_buff_info(buff)
        for ind, value in @kmers:
            backs = get_buffer_level(buff,j,0)
            fronts = get_buffer_level(buff,j,j)

        # clear out buffer state before re-use
        for i in 0 .. len(buff.levels)-1:
            buff.levels[i] = newStringTable() 
            #### thought would be faster, but doesn't work:
            #### clear(buff.levels[i], modeCaseSensitive)

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
        reads_file = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.head"
        (bf,sources,sinks,reals)=load_bf_sources_sinks(reads_file, 10_000)
        bf_cands = get_candidate_paths(reads_file, bf)
    # var read = "ACGTTCGTTTGACACTTCGTTTGTCGTTTGGTTCGTTGTTCGTT"
    # echo reverse(read)
    # echo get_rc(read)
    # echo read # original isn't changed
        
