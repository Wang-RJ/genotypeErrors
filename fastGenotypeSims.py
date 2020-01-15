import numpy, sys

def flag_to_genotype(genoflag):
    geno_dict = {
            0 : '00',
            1 : '01',
            2 : '11'
            }
    return geno_dict[genoflag]

def genotype_to_flag(genotype):
    flag_dict = {
            '00' : 0,
            '01' : 1,
            '11' : 2
            }
    return flag_dict[genotype]

def flag_to_phase(phase):
    phase_dict = {
            0 : 'F',
            1 : 'M'}
    return phase_dict[phase]

# Build all possible biallelic genotype combinations from N individuals
def build_genocombinations(N):
    def base3(x):
        return ((x == 0) and "0") or (base3(x // 3).lstrip("0") + str(x % 3))

    genoflags = [base3(i).zfill(N) for i in range(3**N)]
    genotypes = [''.join([flag_to_genotype(int(x)) for x in geno]) for geno in genoflags]

    return genotypes

# Draws n sites from the same allele frequency based on SFS
class siteGenerator:
    # Upon init, calculate a weight_dict for all 3**3 genotypic combinations for 3 sites
    # Probability weight for each key is sum across SFS
    # Also calculate a normalized SFS (Watterson correction factor) based on sampleSize
    def __init__(self, sampleSize = 50):
        self.sampleSize = sampleSize
        self.weight_dict = {}

        rawSFS = [1/x for x in range(1, self.sampleSize)]
        sumSFS = sum(rawSFS)
        SFS = [freqprob / sumSFS for freqprob in rawSFS]

        def calc_weight(geno1, geno2, geno3):
            def calc_AFprob(geno, i):
                p = i / self.sampleSize
                if geno == '00':
                    return (1 - p)**2
                if geno == '01':
                    return (2 * p * (1 - p))
                if geno == '11':
                    return p**2

            weight = 0

            for i in range(1, self.sampleSize):
                weight += calc_AFprob(geno1, i) * calc_AFprob(geno2, i) * calc_AFprob(geno3, i) * SFS[i - 1]

            return weight

        geno_combos  = build_genocombinations(3)
        exploded_combos = [[geno[i:i+2] for i in range(0,6,2)] for geno in geno_combos]
        geno_weights = [calc_weight(*exploded) for exploded in exploded_combos]

        self.weight_dict = dict(zip(geno_combos, geno_weights))

    # Draw 3 genotypes from N loci
    def draw_genotypes(self, N):
        counts = numpy.random.multinomial(N, pvals = list(self.weight_dict.values()))

        return dict(zip(self.weight_dict.keys(), counts))

# errorGenerator initialized with a matrix of transition probabilities
class errorGenerator:
    def __init__(self, eps01, eps10, eps20, eps02, eps12, eps21):
        self.trans_mat = numpy.array([(1.0 - eps01 - eps02, eps01, eps02),
                                      (eps10, 1.0 - eps10 - eps12, eps12),
                                      (eps20, eps21, 1.0 - eps20 - eps21)], dtype = numpy.float64)

    def emit_errors(self, geno_dict):
        def error_check(genoflag, count):
            stateVector = numpy.zeros(3)
            stateVector[genoflag] = 1

            transitionProbs = numpy.dot(stateVector, self.trans_mat)
            errCounts = numpy.random.multinomial(count, pvals = transitionProbs)

            return errCounts

        # We want to iterate first by each two-character position, not by key
        # This is so that the sampling in error_check is performed once for each genotype call rather than for each genokey
        for call_idx in range(5):
            position_idx = call_idx * 2

            # Emission dictionary empty at beginning of each iteration through position
            emission_dict = {}
            for genokey in geno_dict.keys():
                genotype = genokey[position_idx:position_idx+2]
                genoflag = genotype_to_flag(genotype)

                emittedCounts = error_check(genoflag, geno_dict[genokey])
                for emitflag in range(3):
                    emitkey = genokey[:position_idx] + flag_to_genotype(emitflag) + genokey[position_idx + 2:]

                    if emitkey in emission_dict:
                        emission_dict[emitkey] += emittedCounts[emitflag]
                    elif emittedCounts[emitflag] > 0:
                        emission_dict[emitkey] = emittedCounts[emitflag]

            # Update the genotype dictionary with the emission dictionary after each individual slicing
            geno_dict = emission_dict

        return geno_dict

# build a set of observations by creating a sitegenerator and errorgenerator
# take the observed genotypes and emit errors based on transmission model
def build_related(unrelated_genocounts):
    related_genocounts = {}

    for genokey in unrelated_genocounts.keys():
        P1Geno = genokey[0:2]
        P2Geno = genokey[2:4]
        PartnerGeno = genokey[4:6]

        # Phasedfocals is ordered arbitrarily and counts are drawn for each of the 4 possible inheritance classes
        # Analagous to drawing outcomes from a punnett-square
        phasedfocals = ((P1Geno[0], P2Geno[0]),
                        (P1Geno[0], P2Geno[1]),
                        (P1Geno[1], P2Geno[0]),
                        (P1Geno[1], P2Geno[1]))
        focalcounts = numpy.random.multinomial(unrelated_genocounts[genokey], [1/4.]*4)

        # Nesting additional draws for offspring with each possible inheritance outcome
        # Strategy is to branch outcomes and counts from each genokey

        for phasedfocal,focalcount in zip(phasedfocals, focalcounts):
            phasedoffsprings = ((phasedfocal[0], PartnerGeno[0]),
                                (phasedfocal[0], PartnerGeno[1]),
                                (phasedfocal[1], PartnerGeno[0]),
                                (phasedfocal[1], PartnerGeno[1]))
            phases = (0, 0, 1, 1)
            offspringcounts = numpy.random.multinomial(focalcount, [1/4.]*4)

            # Phasedoffsprings is not ordered arbitrarily because zipped with phase of phasedfocal

            for phasedoffspring,phase,offspringcount in zip(phasedoffsprings, phases, offspringcounts):
                phasedkey = (genokey +
                             flag_to_genotype(sum([int(geno) for geno in phasedfocal])) +
                             flag_to_genotype(sum([int(geno) for geno in phasedoffspring])) +
                             flag_to_phase(phase))

                if phasedkey in related_genocounts:
                    related_genocounts[phasedkey] += offspringcount
                else:
                    related_genocounts[phasedkey] = offspringcount

    return related_genocounts

def writeObservations(out_dict, output_path):
    # Generate all 3**5 possible genotype and phase combinations from 5 individuals
    genocombinations = build_genocombinations(5)
    allKeys = ([entry + 'M' for entry in genocombinations] +
               [entry + 'F' for entry in genocombinations])

    with open(output_path, 'a') as output_file:
        if output_file.tell() == 0:
            output_file.write('\t'.join(allKeys) + '\n')

        for key in allKeys:
            if key in out_dict:
                output_file.write(str(out_dict[key]) + '\t')
            else:
                output_file.write('0\t')

        output_file.write('\n')

def main(output_path, n, segSites, sampleSize,
         eps01, eps10, eps20, eps02, eps12, eps21):

    eps_vec = (eps01, eps10, eps20, eps02, eps12, eps21)

#    print("Epsilons: ", end = '')
#    print(*eps_vec, sep = ', ')

    eps_vec = [float(eps) for eps in eps_vec]

    SG = siteGenerator(int(sampleSize))
    EG = errorGenerator(*eps_vec)

    for i in range(int(float(n))):
        sample_dict = build_related(SG.draw_genotypes(int(float(segSites))))
        out_dict = EG.emit_errors(sample_dict)
        writeObservations(out_dict, output_path)
#        print("Wrote simulation to: " + output_path)

if __name__ == "__main__":
    main(*sys.argv[1:])
