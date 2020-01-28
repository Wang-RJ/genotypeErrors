#!/usr/bin/env python3

import sys, os

def translate_genotype(genotype):
    geno_dict = {
        0 : '00',
        1 : '01',
        2 : '11'}
    return geno_dict[genotype]

def genoCombinations():
    def base3(n):
        return ((n == 0) and "0") or (base3(n // 3).lstrip("0") + str(n % 3))
    for i in range(0, 3**5):
        geno = base3(i).zfill(5)
        yield ''.join([translate_genotype(int(x)) for x in geno]) + 'M'
        yield ''.join([translate_genotype(int(x)) for x in geno]) + 'F'

def push_entry(entry, dictionary):
    if entry in dictionary:
        dictionary[entry] += 1
    else:
        dictionary[entry] = 1

def main(gmother, gfather, offspring, partner, ftwo, minGQ, phase_file, output_path, vcf_path):
    contig_dict = {}
    output_dict = {}
    
    with open(phase_file, "r") as phase_information:
        for line in phase_information:
            l_list = line.rstrip().split("\t")
            contig = l_list[0]

            if contig in contig_dict:
                contig_dict[contig].append(tuple(l_list[1:]))
            else:
                contig_dict[contig] = [tuple(l_list[1:])]

    allKeys = list(genoCombinations())

    with open(vcf_path, "r") as vcf_file:
        header_list = vcf_file.readline().rstrip().split("\t")

        individual_idx = [header_list.index(individual) for individual in [gmother, gfather, offspring, partner, ftwo]]
            
        for line in vcf_file:
            l_list = line.rstrip().split("\t")
            contig = l_list[0]

            if contig not in contig_dict:
                break
            else:
                position = l_list[1]
                
            value_list = contig_dict[contig]
        
            for start,end,phase in value_list:
                if int(start) <= int(position) and int(end) >= int(position):
                    full_info = [l_list[idx] for idx in individual_idx]
                                
                    abbrev_geno = "".join([geno_entry[0:3] for geno_entry in full_info]).replace("/", "") + phase

                    if abbrev_geno not in allKeys:
                        continue
                    else:
                        min_gq_entry = min([int(geno_entry.split(":")[3]) for geno_entry in full_info])

                    if min_gq_entry >= minGQ:
                        push_entry(abbrev_geno, output_dict)
    
    def writeOutput(dictionary, path, keys = allKeys):
        with open(path, "a") as output_file:
            for key in keys:
                if key in dictionary:
                    output_file.write(str(dictionary[key]) + '\t')
                else:
                    output_file.write('0\t')
            output_file.write('\n')
    
    writeOutput(output_dict, output_path)
    
if __name__ == "__main__":
    main(*sys.argv[1:])
