###new scaffold from staples###
### 2019 03 06 ###
### made by Fabian Schneider ###

### script reads all domains and their respective base.sequences to find all skipped and inserted bases. From this new strand domains are created using debruin scaffold sequence and a given set of staples.###

from __future__ import (absolute_import, division, print_function, unicode_literals)

import os
import sys


#copied from nanodesing/strand-statistics.py
# The following imports are designed to try and get the nanodesign package
# imported regardless of whether or not you have it installed in the
# site-packages. If you have it in site packages, the first block should just work.
# Otherwise, it will assume you can import it based on a relative path from this source
# file's directory, and try to do so by adjusting the system paths temporarily.

try:
    import nanodesign
except ImportError:
    base_path = os.path.abspath( os.path.join(os.path.dirname(os.path.abspath( __file__)), '../nanodesign/'))
    sys.path.append(base_path)
    import nanodesign
    # If the import fails now, we let the exception go all the way up to halt execution.
    sys.path = sys.path[:-1]

from nanodesign.converters import Converter
from nanodesign.converters.dna_sequence_data import dna_sequence_data


def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter


def rev_comp(staple_sequence):
    """ Generate reverse complement """

def create_staple_sequence(staple_length):
    """ Generate new staple from scaffold sequence """

def get_strand_length(strand): # computes total length of strand in bases
    length = 0
    for base in strand.tour:
        length = length + 1 + base.num_deletions
    return length

def rev_complement( seq ):
    comp = {"G":"C","C":"G",
            "A":"T","T":"A",
            "g":"c","c":"g",
            "a":"t","t":"a",
            "N":"N","n":"n"}
    rev = []
    for i in xrange(len(seq)):
        rev.append(comp[seq[len(seq)-i-1]])
    return "".join( rev )


def generate_staple_from_debruin(debruin_sequence, staple_length, scaffold_start):
    new_staple_sequence = debruin_sequence[scaffold_start:(scaffold_start+staple_length)]
    scaffold_start+=staple_length
    return new_staple_sequence

def pick_provided_staple(scaffold_reservoir, staples, staple_length, dummy_position):
    found_staple = False
    for position_index, item in enumerate(staples):
        if len(item) == staple_length:
            return staples.pop(position_index)
            found_staple = True
            break
    if found_staple == False:
        return generate_staple_from_debruin(scaffold_reservoir, staple_length, dummy_position) #pCS2
        dummy_position += staple_length




#pCS2 = "ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGGTAGAATTGGTAAAGAGAGTCGTGTAAAATATCGAGTTCGCACATCTTGTTGTCTGATTATTGATTTTTGGCGAAACCATTTGATCATATGACAAGATGTGTATCTACCTTAACTTAATGATTTTGATAAAAATCATTATAATGCCAGGGTCAGTGAGATGTAACACTACAGACCTCAGCTTCCTCTACTTTGCCTGGCATAGCCTACTGCACCATCCTTTTCAATGGAAGTAATTTGTCCCCCATGGGTTACACCCCCACCTAGACTACTTAAGAAGGTGCTATTGCTCCTGATTTCCAATGTTTGAATCTCCCTCCCACCATACCAGGCATAAAGGCTGTCACTCCCCAATACTAAGTGGCCCCAAGTTCACCACAATACAGCCAAGTGCATTGTTATCACATTCAGTCAAAGAATTTTCATTGGCAACAGTCCTCTCCAATCTGCCTATTCACTGTCCCTACATCAGGTGATAGTCCAGTTGTCAGACATGTGATCCCTGTTTGGGGACCCACAAGGCCTATCCAACCCATTGCTTTTAAAAGAGCTGGGGTTCTGGTGACTTCCACAACAAGCCCATCTGCTCTTAGTACAGCTTGTCTTAAAATTAGCTAGACACTGTATCATACCCCCTCCACTGCCCTATGATAAGGACAGCAGAGGGTTAGATTTACTTATCTGCACTCATGTGCCCCCAGCCTGGTACCAAAACTGATAATCTACTCCACAGCTGCATACCTGCCATCTTAGCCAAACCACACCAATCCAATTGCTGGGTCTCATATCATCTCTGAAGAAGTGATTACCTCATTACATGCTTGCTTGATATACCCACTTTCCTTATGAGGCACACCCATACTGAGGTAGCTATTTTACTATTGACTGGGTAAGCCAGTACAAAATAATCATAGGTGCACAGTCTGGACATGGGAACCCCTCTTCCTACTCACCTCTCTCAGAGGCAATGCCCCTTCATTAACCATGCACCCTAGGCCTTGGTACAAGTTGCCAACCTGCAGAGAACCAGTTAGCTGAAAGGCCAGGTTCTTCTTGATCTTGGTGCTTAGGGCAAACTTTTGGGATAGGCAAAAACAAGTCCTTTACACACCTCCAGTCTACAGTTCCATCAAAGTGGAACTTGGAATGTTATAAACTGCATCTATACTGTACTGCTAGTAGGCTTGTTTAATGCTTTCAACTCCTAGAGAATCACTTCAATCCTACCAACCAATAGTTCCCCATCACACAGATGGCTAAAGACCCTACCTCCCCTAAAGCCCTACTTCACACTGCAGGGACCTCCTAAGGGTCTTGGACCCTCCTCACTTGTTACTAGGCTATCATTCTACTGGTGCAGATACAGATTTTAGGTACCTGTCCTTGATGCACTTAGCTTTAACATGTCACATGCCAATGAAAATAGCAGCCACAGGGGGATATTCTAGGAGCAACATGGAGACTGTCTGCTACAACCACCAAACAGCTCATCCTCTGCATGAATACTTACTCTAGACCCCTGCATTCACCTTCCATGATTAATCACCAACTGAGTCTTGCAAAACCCAATTATGGTTATAGACAGACTGGCCCATGTTACCACTTAAACTCCATCTAGTCTGCAGTGTTCTGACCTATGCCTGAGGGTCCAGATAACCACTCTTCACTAAGCCCCCTGGAAGACAGTTAAGTCTATCTATCAGCCTCACATCCATCCCATTCTCTCTTTGAGATCTTATCAGGACTGATTGTATGTGAACAATATACTCCCAGGCTGGCCATTAGAGATTGGTTAGTGGTTGGCTGGACCAACACCCTTAAGCAACTAATAAGACCATGTGGGGTGTCCAGGCCTCATAAATCCTGAGCCTCTCAATCTTCCAGTGGCTCACCATTCCCATCCACACAACCTATTGGAACCATCAGAGTCATCACCTATATAGGGTGAAGTAGAGCTCCACCACCCCATAACTATACAGTGCACCTGGCCTCTGAGACTTGCTGTGATAACACAAACCCCAGGGCCAAATTGGGAAGATGCAATGACATTCCATTACTTCTCCAGGACAACATTTATTTGTTGCTGACTCCAGACTTAACAAGGAGTATTTGAGGCCACTGTTGGGGCCCTAACCCACCCAGCTCTTGACAGCTAATCTTAACCTCAAGGATTGACCATTGGGGGAGTCAACCAGCTTTGGCCCTGCTCAAGCAAACAATCCCCAGATGACCTGCTGGCAGGCACCCCTAGAAATATTGTGGGTTGCACTATTTCTTGCTCAGATGCTTCCCACTATAGCCCCACTCCAACTAGTACCACATCATTTTATGGATCAATCACAAAGGAGCCCCTATAGGCTGATTCTGTGTTTCAAAAAGTGAAGCCTATGTGCTCCATTGAGTCAGTAAATTAACACCACTAGAAGTTATCCAGCTAGTCAAGGGTTTTCCATATGTCCCAAGGTACTAGCTACAGCACTTGAATATGGTCTACTATGAGTGGGACTATTAATGTCATCTTTCTGTATAGCTGCCCCATTTAAAGCATATATGAAGTGCTAACATTGTGTTGGATCTTTATGCCACCAGATCCCCCTTACAGTATTCTCACCCAAATGGGATGAGAGTTGTATTAAGGCTAGTTAATCCCTAAGTCCATACATAGGACCACAGTAGCTTGGCTATGCTCACAATTTGCCATGTATATCAGTCCCAGAGGTATGCAAGACATTGCCCACATTTCAGAATGGTCCTGTCTTGTGGTTCCTTGCCTTTTGCACATTACCCTTGGCCTTCTACATAATGTGAGTCCTGGCTGCAATCTATGTATTGAAAGTAGGTAACCCTTCTCATTGTCAAGAACATCTGTTCCTAGCCACCCTGACCCTGCCTCCCTGAAGCTGCTCCCATAGGCCCACTGATCTCACTGCTTCTAATTGCATCCCTTGACTCTCTATCCCCTGAGAATTCCTGCTTGTGCAGCCTTGTGAGCCACTACCCAGTCAGGGTGTTAAGCCTTACCATGACCAGGTGCCCTCTATGACTAGCCCAAGACTCACAGACAGGATAGCCAGCTGTACCCCACACTCCTCTTGCCCTGGGCACTCTCATCAGTACCCTCAATAACCTAGCTGTCCACTCACTATCCTCAGACCAATTCTCCTTCTGCCAAAGGTCCAACATAGAGGACTCATCTGGGCTCAGGCAACTCAGGGCATATGGGACCATAGCTCAACTACAATAAAACCAGGGGCCACATAAGCTAGCAATCAGTGGGTGTAAGTGTCTGGGTGCATCAGCTGGTAAGACTGCCAGATTCTATGCAGGCCAATCATGAACCTACTAGTCCCTCATGGTCACTTTGTTCATCTAACAGATAGTTGCTTCAAACTCAAGACCTTACTAACCTTAGTCACTGGTTCTACCAGCCCAGGAATCCTCCTTGGGCTAGCCTCCAAGAAGCTTCAGACTCCTTAATACCAAGCAGCAGGAACAGTTTCAGCACCAAGAGGTGAGAACTGCTGTTAGCATCACTGACTTACAATCTCATGATAGCTAAGCATTCTTAGAACTCTTATATAACTCCCTATTACCAGAGACAGCCCTCAGGATTACTGATGGACCTTCTTCAAGTCTTTTCTAGCACAGGACCTACAGGGAGAGGGAAAGTCAGCTCCTATGGGGAATACAAGAGTTAGTCTTCAGTTTGTGGCTGAGGCTTCTTTGTGCTGGTTGCCTCTAAGACAAATTATTCTGGAATCTAATCAAGTTTACAAGGTCAGCATTATGACCCAGGTCAATAGCTTCTGGGAGGCTCTTCTAGTTTCCACCTTTGTACCAGTATCCTATTTAAGATAGGTCATGCATAGGGCTCTACACAATGAGTTCAGAGCCTTTGCAGCTTAAGGAAATCTCTCCTCCATAAAAGTTCTGCTTTGATGTCAATTCCCCTTGTACTTCCTGTGACTGCAACACAGTTGGACTGGTATTCCTCAACATCAACAATTAGTAAAAGCAATGGTTGTGACCAAGGAATGAGACCAGCAATTAATTGACATAGTTTAACTGTCAAACACTGAGAGATGGAGCCTGCAAGTGGGGCTTATCCCACAGAGTACTGTTTTCTGAGGAGCTGAGTTAACAGGTGGAGAAAATCCCAAAGATACTCTCCCAATGGCCAGTCCAAGCTTGAAGTCAATGCAACCTCTAGTAATAGGTAGAGGCTGCCACAAATGATGCCTACAAAGAGTAAACTTCTATATGTGGTAATGGCATCCAAGTAAAGGGCTAATGCAGAACAGAGCTTTTGTTAATGGGCAAGCTACTACTAAACCCTCTGGAGCATGTGTTAGTAGATACCTTCAACCTTGGAGCTCACTCAGAAAAAGCCTGAACCCAGAATTAAAACATCCTGGTGTACTCTGGCACCAGTGACACCATGGCTTTTTCAGGCCCCTGTGCCTGTTGCATGCAAACCTAAGAGTCTAAAGGTAATCTGTCAGCCAGGAGAGTGGTCTCCATGCCTCAGTGTGTAGGTGGGCAGTCATTCCTATACCATATAAAATGGAGTACATCTACAAGCACTGGACTACCTAGGTTTCCTGGATAGTAACTCTAACCAGACACACATGGCCTGGAGTGGAGTTACCCAACAGAACCCTGTCATATTTCCTAATCCACTAAAGTCTGTTACATACAGGCAGTTTTTGCTCTCTGCTGCCTGCCCAGCAGCTATGGTGATCAACTGGCTAGAATGACTTGTGTACCTACCCTAATTTTAATTTCACTCTGACATGCAGTTCAACAGGCTACCCCTTTTTATTGGTGGCTACATTTGGCTTAATTCTAACTGGTCAGAGATAAAGCTCAGCAGTGATATGAGCCCTTAGGAAAATTCCAGAGCACACACTAGCAGATGTTTTATCATGGAATATTATCAATGTAGCTCTGGGGGTAACTTTGAAAAGCTGTTCACTTACCTTTCCAGGGATAACTTGCATTTAGACATACTAGATTATCTTACATTCTGATATCACTAGTGTCCTGAAACAGGGCTTCATCATAATCCTTGTCACCTGTTAAACCTTTTATATGCCCTTGCTAGGGCCCAACTCTGTCTCCTGTAAAGTTAAAGTGCCTTCAGGAAGTTGGTGAGTAGACCTGAATAATTCATACTTTATCTCCTACACCAGGATGACTCAGTAGGACATCTCAGTCTCTAATACTCATACAAATACACTCTATTAGGGGCTGGGACACCTGATCCAGAAACAAAGTTTCTTAAGTACCTCTGTACAATTGAAGATAATGGTACTCAATGATATTTTGTGATTCAAATTTACAGAGAGCAGTAATGTTCCCAGTTTATCCTTAGATAGAGCAAAAGATTAAAGAACCTGGTTAATATTTAGGATCTGATGATCCTGCACACTTCTGAATGCCTTGAAATTCTTGTCCATTCAATACCCTGGTCTTAGGCCATATCCTGTTCTAGAACACTTATGCACAATCAATATGCATGGGCTGTGTGAGGGCTGCTATCTGTAGGCACTAATGGATTTGTAAGCTGGCTCTGTTTCTCTTACTTTCAGTATACCTATCAAGGTTTTGGTGTTGAAACCTGACTGACCACTGGCAAAGGGTGGGAAACTTAGGTTCCCTAGTATGCTAGCTTAGCAACCATAAGGTCTATTGTCTCAAATCAAAATCTAGGGGAACATGATCTGGATGCCCATATTCAAGCCAACAAATCTTGAACTGGGATTCCAAGGGCCTAACAACTATTCCACTTGGCACTTTAGCAGGCTCATTCATGCTCTGATCAGTTGCAATTTATGAAACTGTAGTCTCAGGTCCTAGGAATTATAGGTTGTACACTGGGCCTCCTGGGAATTGTTTGCTACTTGACCTCTTTCATGTACATGTTGATGGGGTAGTACTAATTAAGAGCCATAGTCATAAGTAGTTGGGCATGCTGAATCATCCAGTAGAAGGCTCCAGCATGAGGTGGCATTGATGAACAGCCTAAGCTTATTGCAAGCATGGCACATAGCACCTAACTTCATGACTGTTATTCATTTACCTGAGTATAGGATTTCATATAGATCCAAATAGGAACTAAGGTTCATAGTAGCCATCATGTTCAAAGCTATAACAAACATTAAATATATCCATGTCCTAAATGCATTAGGTGAACCAAGTCACAGGTACAGAAAGGGGATGCAGCACAAAAATTGATACTTGTAGAAAGAGGATTCAGATATATTGAGCTGTGCTACCATTAATAAAGATGGGTCCCATGAAGAGCATCTTCTGTTGTTAGGTCTGACTATGGCAGCAAATCCAGGTAAGGCAAGGACTAGTTCTCAAGTGACCTTTAAGTTGACTTTTAGCCTGTATTTCAAGGCTTAGAGGGCACAACTGTGTCTATAAATTGCCTAATGAATCCATAGACTTTACCAATGCTGTCTAGATATGCTTAAAGGAAGAGTGTTGTAGTGTATGAACTTTCTATTTGCTGCTGAGATATTAGTTGTGGATCCTAGTTGAATTTCTAAGTTAGGCAGGTATCTGGCTTGCAGGATCCATTTCTGCAAAGTAAGAGGGGGTCATTGCAGACCCATGCTATGTCTGTGAAAGATCTGTGCATATCTAGAGGTTAGGACTTCAGCTACCTGGGTTTGCATAAGAATAAGCAGACTAAAAGGCATGACAGAAGACTTCTTATGGGTAGCAGTCTTATTTTTAGAGTTCTTGGGAGATCAGATCTAAGGAGATTCACAAGAAACTAGAGTAGCATTGGAGATGAGTAACATATACACATATTACTACATGGTTTACTGCCTTAAATGTGGCCTAGGGTATTATGCTGGAGAGCTAGGTGTCAGTTCTAAAAAGGATAAGTCATGGGGGCCTTTCACAGCAAGTTACAAACTATGTTTAGTATCAGAACTTATTAAACATAACAGTAACCAAAGCCATGGTATGATGTACTATCTTTTTGTAGCACTAGGATACTATATTAACTACTCTTTAGTTTTGTCATGAGAAGTACACCTTATAGCATAGATTGAATGGGTGGATGTCTCTGTGGAGGTGTGTGCTTGGATTAGATGGTGTGAATAGAAAACTAAATACTGGATTGGACAAGTGTTTGTCTAATGTATGGACAGTGGATAAAC"


def main():


    # create path to caDNAno file
    current_working_directory = os.getcwd()
    file_name = input("Cadnano design file to load (format: \"test.json\"): ")#sys.argv[1]
    file_name = os.path.normpath(current_working_directory + "/" + file_name)
    print('file_name to be opened', file_name)


    # set sequence to assign to scaffoldself.
    scaffold_seq_name = input("Scaffold as template (format: \"pCS2\", set as standard): ")#sys.argv[2]


    #read provided staple strands
    file_name2 = input("Textfile with staple strand sequences (format: \"providedstaples.txt\"): ")#sys.argv[3]
    file_name2 = os.path.normpath(current_working_directory + "/" + file_name2)
    print("file to be opened ", file_name2)
    with open(file_name2) as f:
        provided_staples = f.read().splitlines()
    with open(file_name2) as f:
        #second copy of provided_staples for later use
        provided_staples2 = f.read().splitlines()
    print(provided_staples)






    # read cadnano file and create dna structure.
    converter = read_file(file_name, scaffold_seq_name)
    dna_structure = converter.dna_structure

    # determine domain information
    dna_structure.get_domains()


    # determine scaffold strand id, alert if multiple scaffold strands
    # prints starting location if not circular strands

    # strand id of last found scaffold strand
    scaffold_id = -1

    for strand in dna_structure.strands:
        if strand.is_scaffold:
            if scaffold_id != -1:
                print('multiple scaffold strands!')
            scaffold_id = strand.id
            if not strand.is_circular:
                print('scaffold strand starts in helix ', strand.tour[0].h, ' base number ', strand.tour[0].p)
            else:
                print('scaffold strand is circular, passes through helix ', strand.tour[0].h, ' base number ', strand.tour[0].p)






    # test to replace staple sequences
    scaffold_reservoir = dna_sequence_data.get(scaffold_seq_name, None)
    dummy_sequence = scaffold_reservoir
    domain_length=0
    dummy_position=0
    staple_pos=0
    base_step=0
    final_sequence=""
    for strand in dna_structure.strands:
        base_step = 0
        if strand.is_scaffold:
            for domain in strand.domain_list:
                loop_bases = 0
                domain_length=len(domain.sequence)
                for base in domain.base_list:
                    loop_bases += base.num_insertions
                if domain.base_list[0].across==None:
                    skipped_bases = domain.sequence.count("N")
                else:
                    skipped_bases = 0
                domain_length -= skipped_bases
                domain_length += loop_bases

                if domain.base_list[0].across==None:#.seq=="N":
                    dummy_sequence = ""
                    dummy_sequence = generate_staple_from_debruin(scaffold_reservoir, domain_length, dummy_position)
                    for i in range(skipped_bases): # doesnt work if you skip all bases of a segment!!!
                        domain.base_list.pop(-2)
                    # this should implement looped bases to the base_list
                    for i in range(domain_length  - len(domain.base_list)):
                        domain.base_list.append(domain.base_list[-1])

                    for n, item in enumerate(dummy_sequence):
                        domain.base_list[n].seq = dummy_sequence[n]

                    final_sequence=dummy_sequence
                    domain.sequence = final_sequence




        if not strand.is_scaffold:
            staple_length = get_strand_length(strand)
            staple_length = 0
            #recalculate staple length including insertions and deletions
            for domain in strand.domain_list:
                loop_bases = 0
                domain_length = 0
                staple_length += len(domain.base_list)
                domain_length = len(domain.base_list)
                double_stranded = True
                bad_base_indices = []
                for i, base in enumerate(domain.base_list):
                    if base.seq == "N":
                        if not base.across == None:
                            staple_length -= 1
                            domain_length -= 1
                            bad_base_indices.append(i)


                    if base.num_insertions != 0:
                        loop_bases += base.num_insertions
                        staple_length += base.num_insertions
                        domain_length += base.num_insertions
                        if not base.across == None:
                            for i in range(base.num_insertions):
                                dna_structure.domain_list[domain.base_list[0].across.domain].base_list.append(dna_structure.domain_list[domain.base_list[0].across.domain].base_list[-1])
                                dna_structure.domain_list[domain.base_list[0].across.domain].base_list[-1].num_insertions = 0
                                domain.base_list.append(domain.base_list[-1])
                                domain.base_list[-1].num_insertions = 0
                                domain.sequence = "".join(base.seq for base in domain.base_list)
                                dna_structure.domain_list[domain.base_list[0].across.domain].sequence = "".join(base.seq for base in dna_structure.domain_list[domain.base_list[0].across.domain].base_list)

                        else:
                            for i in range(base.num_insertions):
                                domain.base_list.append(domain.base_list[-1])
                                domain.base_list[-1].num_insertions = 0

                    if base.across == None:
                        double_stranded = False

                    bad_base_indices = [x+loop_bases for x in bad_base_indices]
                for i in reversed(bad_base_indices):
                    domain.base_list.pop(i)
                    final_sequence="".join(base.seq for base in domain.base_list)
                domain.sequence = final_sequence

                bad_base_indices = [domain_length+len(bad_base_indices)-x-1 for x in bad_base_indices]
                for i in bad_base_indices:
                    final_sequence="".join(base.seq for base in dna_structure.domain_list[domain.base_list[0].across.domain].base_list)
                    dna_structure.domain_list[domain.base_list[0].across.domain].base_list.pop(i)
                    final_sequence="".join(base.seq for base in dna_structure.domain_list[domain.base_list[0].across.domain].base_list)


                if not base.across == None:
                    dna_structure.domain_list[domain.base_list[0].across.domain].sequence = final_sequence
                #if double_stranded == True:
                    #print("the two domains: ", domain.sequence, dna_structure.domain_list[domain.base_list[0].across.domain].sequence)




            dummy_sequence = pick_provided_staple(scaffold_reservoir, provided_staples, staple_length, dummy_position)
            dummy_position_local = 0
            for domain in strand.domain_list:
                    domain_length = len(domain.base_list)
                    if domain.base_list[0].across != None:
                        skipped_bases = domain.sequence.count("N")
                    else:
                        skipped_bases = 0
                    domain_length -= skipped_bases
                    domain.sequence=dummy_sequence[dummy_position_local: domain_length+dummy_position_local]
                    dummy_position += domain_length
                    dummy_position_local += domain_length
                    for i, base in enumerate(domain.base_list):
                        if not base.seq == "N":
                            base.seq=dummy_sequence[base_step]
                            base_step+=1
                            if not base.across==None:
                                dna_structure.domain_list[domain.base_list[0].across.domain].base_list[-i].seq = rev_complement(base.seq)


                        else:
                            base.seq=dummy_sequence[base_step]
                            base_step+=1
                    if not domain.base_list[0].across == None:
                        dna_structure.domain_list[domain.base_list[0].across.domain].sequence = rev_complement(domain.sequence)

        loop_bases = 0
        domain_length = 0



    output=open("created_scaffold.txt", 'w')
    final_sequence="".join(domain.sequence for domain in dna_structure.strands[0].domain_list)
    print("created new scaffold sequence\n\nScaffold sequence:", final_sequence, file=output)
    output=open("created_scaffold.txt", 'a')
    print("\nstaples used:", file=output)

    used_staples = []
    created_staples = []
    for i, strand in enumerate(dna_structure.strands):
        final_sequence = str("".join(domain.sequence.strip() for domain in strand.domain_list))
        if final_sequence in provided_staples2:
            used_staples.append(final_sequence)
        else:
            created_staples.append(final_sequence)

    print("\nused staples from provided list:\n", file=output)
    for strand in used_staples:
        print(strand, file=output)
    print("\nnew created staples with ", scaffold_seq_name, " :\n", file=output)
    for strand in created_staples:
        print(strand, file=output)

if __name__ == '__main__':
    main()
