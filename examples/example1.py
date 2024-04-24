#!/usr/bin/env python
"""
Example usage of HGVS library.

To run the script, first begin by opening a terminal in the root directory
of this software package. Note, the root directory should contain
`setup.py`.

Second, obtain genome sequence in FASTA format, which is required in
example. Genome sequence can be fetched using the following commands:

  cd /tmp
  curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
  tar zxvf chromFa.tar.gz
  cat chr*.fa > hg19.fa
  rm chr*.fa chromFa.tar.gz

This example script can be run using:

  python examples/example1.py

The following output should be displayed:

  chr11 17496508 T C
  NM_000352.3(ABCC8):c.215A>G
  ('NM_000352.3', 'c', '>', CDNACoord(215, -10), CDNACoord(215, -10), 'A', 'G')

"""
from __future__ import print_function

from __future__ import unicode_literals
import pyhgvsv as hgvsv
import pyhgvsv.utils as hgvsv_utils
from pyfaidx import Fasta

# Read genome sequence using pyfaidx.
genome = Fasta('/tmp/hg38.fa')

# Read RefSeq transcripts into a python dict.
with open('pyhgvsv/data/genes.refGene') as infile:
    transcripts = hgvsv_utils.read_transcripts(infile)


# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)


# Parse the HGVS name into genomic coordinates and alleles.
chrom, offset, ref, alt = hgvsv.parse_hgvs_name(
    'NM_000352.3:c.215A>G', genome, get_transcript=get_transcript)
print(chrom, offset, ref, alt)
# Returns variant in VCF style: ('chr11', 17496508, 'T', 'C')
# Notice that since the transcript is on the negative strand, the alleles
# are reverse complemented during conversion.


# Format an HGVS name.
chrom, offset, ref, alt = ('chr11', 17496508, 'T', 'C')
transcript = get_transcript('NM_000352.3')
hgvs_name = hgvsv.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript)
print(hgvs_name)
# Returns 'NM_000352.3(ABCC8):c.215A>G'

# Format an HGVS name for a structural variant (deletion).
chrom, offset, ref, alt, sv_length = ('chrY', 24861625, '', '', -4780)
transcript = get_transcript('NM_001388484.1')
hgvs_name = hgvsv.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript, sv_length)
# Returns 'NM_001388484.1(DAZ4):c.1210-436_1354-437del4780'

# Format an HGVS name for a structural variant (insertion).
chrom, offset, ref, alt, sv_length = ('chr17', 8141778, '', 'TTCTCCCCCCTTGAACTTGAGCTCAATTC', 29)
transcript = get_transcript('NM_002616.3')
hgvs_name = hgvsv.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript, sv_length)
# Returns 'NM_002616.3(PER1):c.3600+26_3600+27ins29'


hgvs_name = hgvsv.HGVSName('NM_000352.3:c.215-10A>G')
# fields of the HGVS name are available as attributes:
#
# hgvs_name.transcript = 'NM_000352.3'
# hgvs_name.kind = 'c'
# hgvs_name.mutation_type = '>'
# hgvs_name.cdna_start = hgvs.CDNACoord(215, -10)
# hgvs_name.cdna_end = hgvs.CDNACoord(215, -10)
# hgvs_name.ref_allele = 'A'
# hgvs_name.alt_allele = 'G'

print((hgvs_name.transcript,
       hgvs_name.kind,
       hgvs_name.mutation_type,
       hgvs_name.cdna_start,
       hgvs_name.cdna_end,
       hgvs_name.ref_allele,
       hgvs_name.alt_allele))
