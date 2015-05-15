# not explained well in docs:
# 1 - how to return array of strings (e.g., get k-mers)
# 2 - access global table/array in function (used in get_rc)

import tables
import bloom
import strtabs

const
    bases = ['A', 'C', 'G', 'T']
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}.toTable()
    k = 27
    fp = 0.01
    j = 10
    read_len = 100

proc get_kmers(r: string, k: int, kmers: var openarray[string] ) =
    var i = 0
    while i < len(r)-k+1:
        kmers[i] = r[i .. i+k-1]
        i = i+1

proc reverse(s: var string): string =
    result = newString(len(s))
    for i in 0 .. s.high:
        result[i] = s[s.high - i]

proc get_rc(dna: var string): string = 
    result = reverse(dna) 
    for i in 0 .. dna.high:
        result[i] = complements[result[i]]

# proc get_canons()
   
proc load_bf_sources_sinks(fname: string, numreads: int) =
    discard """ loads k-mers into bloom filter
        loads to lists ends of reads as potential
        sources and sinks (or nodes j away from sinks)
    """
    var sources = newStringTable() #newSeq[string]()
    var sinks = newStringTable() #newSeq[string]()
    var reals = newStringTable() #newSeq[string]()
    var kmers: array[0..read_len-k+1, string]

    var bf = initialize_bloom_filter(capacity = numreads, error_rate = fp)
    var f_hand = open(fname)
    var line_no = 0
    for line in f_hand.lines:
        if (line_no + 1) mod 10_000==0:
            echo ($(line_no+1) & " read k-mers processed " & 
                $(len(sources)) & " " & $len(sinks) & " " & $len(reals))
        get_kmers(line,k,kmers)
        sources[kmers[0]] = nil
        sinks[kmers[read_len-k]] = nil
        for i,value in @kmers:
            if value!=nil:
                # echo i, " ", value
                bf.insert(value)
                # reals only for debugging:
                reals[value]=nil
        inc(line_no)

        # continue
        # process(line)
    f_hand.close()


# def load_bf_sources_sinks(filename,j,numreads):
#     """ loads k-mers into bloom filter
#         loads to lists ends of reads as potential
#         sources and sinks (or nodes j away from sinks)
#     """
#     sources = []
#     j_sinks = []
#     reals = []
#     B = BloomFilter(capacity = numreads * (read_len-k), error_rate=fp)
#     line_no = 0
#     with open(filename) as f:
#         for line in f:
#             if line_no+1%100000==0:
#                 print str(line_no+1) + " read k-mers processed"
#             read = line.rstrip()
#             sources.append(read[:k])
#             sources.append(get_rc(read[-k:]))
#             j_sinks.append(read[-(k+j):-j])
#             j_sinks.append(get_rc(read[j:j+k]))
#             kmers = get_kmers(read,k)
#             # insert canonical (lex-min) kmers only
#             canons = get_canons(kmers) 
#             B.update(canons)
#             reals.extend(canons)
#             line_no +=1
#     reals = set(reals)
#     sources = set(get_canons(sources))
#     j_sinks = set(get_canons(j_sinks))
#     return (B,sources,j_sinks,reals)


when isMainModule:
    var reads_file = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.100k"
    load_bf_sources_sinks(reads_file, 100_000)
    # var read = "ACGTTCGTTTGACACTTCGTTTGTCGTTTGGTTCGTTGTTCGTT"
    # echo reverse(read)
    # echo get_rc(read)
    # echo read # original isn't changed
        

    

# # if ACGT line get k-mers as list
# # store first, last, last - 10 k-mers 
# # put all in BF
# # import bloom
# echo(bf)                                               # Get characteristics of the Bloom filter
# echo(bf.lookup("An element not in the Bloom filter"))  # Prints 'false'
# bf.insert("Here we go...")
# assert(bf.lookup("Here we go..."))


# def get_kmers(r,k):
#     return [r[i:i+k] for i in range(len(r)-k+1)]

# def get_rc(dna):
#     rev = reversed(dna)
#     return "".join([complements[i] for i in rev])
    
# def get_canons(kmers):
#     return [min(x,get_rc(x)) for x in kmers]