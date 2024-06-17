import argparse
import pandas as pd
from Bio import SeqIO
import logging
import re
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, NotI, SacI
from Bio.SeqRecord import SeqRecord

def find_transcript_sequence(plasmid_sequence, promoter, terminator):
    logging.info("Finding transcript sequence.")

    """Extract the transcript sequence from the plasmid based on promoter and terminator in a circular DNA context."""
    plasmid_sequence = plasmid_sequence.upper()
    promoter = promoter.upper()
    terminator = terminator.upper()
    
    circular_sequence = plasmid_sequence + plasmid_sequence
    
    start_index = circular_sequence.find(promoter)
    if start_index == -1:
        raise ValueError("Promoter not found in the plasmid sequence.")
    
    transcript_start_index = start_index + len(promoter)
    
    polyA_index = circular_sequence.find(terminator, transcript_start_index)
    if polyA_index == -1:
        raise ValueError("Terminator not found in the plasmid sequence.")
    
    transcript_end_index = circular_sequence.find("CA", polyA_index + len(terminator))
    if transcript_end_index == -1:
        raise ValueError("End 'CA' of the transcript not found after the terminator.")
    
    transcript_start_index %= len(plasmid_sequence)
    transcript_end_index %= len(plasmid_sequence)

    # circularity logic
    if transcript_end_index > transcript_start_index:
        return plasmid_sequence[transcript_start_index:transcript_end_index + 2]
    else:
        return plasmid_sequence[transcript_start_index:] + plasmid_sequence[:transcript_end_index + 2]

def simulate_ligation_cloning(transcript, insert, analysis):
    logging.info("Simulating ligation cloning.")

    sacI_sites = analysis.with_sites()[SacI]
    notI_sites = analysis.with_sites()[NotI]

    if not sacI_sites or not notI_sites:
        raise ValueError("Required restriction sites not found in the transcript sequence.")

    if len(sacI_sites) > 1 or len(notI_sites) > 1:
        raise ValueError("More than one restriction site found for SacI or NotI.")

    sacI_site = sacI_sites[0]
    notI_site = notI_sites[0]
    # print(f"SacI site: {sacI_site}")
    # print(f"NotI site: {notI_site}")
    if sacI_site > notI_site:
        sacI_site, notI_site = notI_site, sacI_site

    # break transcript apart at restriction sites
    part1 = transcript[:sacI_site]
    part3 = transcript[notI_site:]


    modified_insert = insert #[len(SacI.site):-len(NotI.site)]
    ligated_sequence = part1 + modified_insert + part3

    return ligated_sequence

def get_transcript_parts(transcript_sequence):
    logging.info("Getting transcript parts.")

    """Extract UTR5, CDS, and UTR3 from the transcript sequence."""
    transcript_sequence = transcript_sequence.upper()
    
    # Find the start of the CDS (first ATG)
    # this ignores the idea of uAUGs -- although the length checker will likely raise an error
    cds_start = transcript_sequence.find("ATG")
    if cds_start == -1:
        raise ValueError("Start codon 'ATG' not found in the transcript sequence.")
    
    # Find the end of the CDS (first stop codon: TAA, TGA, TAG) in the correct reading frame
    stop_codons = ["TAA", "TGA", "TAG"]
    cds_end = -1
    for i in range(cds_start + 3, len(transcript_sequence), 3):
        codon = transcript_sequence[i:i+3]
        if codon in stop_codons:
            cds_end = i + 3
            break
    
    if cds_end == -1:
        raise ValueError("No stop codon found in the transcript sequence.")
    
    # Extract UTR5, CDS, and UTR3
    UTR5 = transcript_sequence[:cds_start]
    CDS = transcript_sequence[cds_start:cds_end]
    UTR3 = transcript_sequence[cds_end:]
    
    if len(CDS) < 200:
        raise ValueError("CDS is shorter than 200 bp.")
    
    return UTR5, CDS, UTR3


def main(plasmid_gb, fasta_file, output_file, add_sites):
    logging.info("Starting main function.")

    PGK_PROMOTER = "ctaccgggtaggggaggcgcttttcccaaggcagtctggagcatgcgctttagcagccccgctgggcacttggcgctacacaagtggcctctggcctcgcacacattccacatccaccggtaggcgccaaccggctccgttctttggtggccccttcgcgccaccttctactcctcccctagtcaggaagttcccccccgccccgcagctcgcgtcgtgcaggacgtgacaaatggaagtagcacgtctcactagtctcgtgcagatggacagcaccgctgagcaatggaagcgggtaggcctttggggcagcggccaatagcagctttgctccttcgctttctgggctcagaggctgggaaggggtgggtccgggggcgggctcaggggcgggctcaggggcggggcgggcgcccgaaggtcctccggaggcccggcattctgcacgcttcaaaagcgcacgtctgccgcgctgttctcctcttcctcatctccgggcctttcgacctgcagccc"
    SV40_POLYA_SIGNAL = "cagacatgataagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggagatgtgggaggtttttttaagcaagtaaaacctctacaaatgtggtaaaatcgaattttaaca"

    plasmid = SeqIO.read(plasmid_gb, "genbank")
    plasmid_sequence = str(plasmid.seq)
    oligo_seqs = SeqIO.parse(fasta_file, "fasta")

    transcript_sequence = find_transcript_sequence(plasmid_sequence=plasmid_sequence,
                                                   promoter=PGK_PROMOTER, terminator=SV40_POLYA_SIGNAL)
    transcript_sequence = Seq(transcript_sequence)

    inserts = []
    transcript_names = []
    for oligo in oligo_seqs:
        if add_sites:
            NotI_site = NotI.site
            SacI_site = SacI.site
            insert = SacI_site + str(oligo.seq) + NotI_site
        else:
            insert = str(oligo.seq)
        inserts.append(insert)
        transcript_names.append(oligo.id)

    rb = RestrictionBatch([NotI, SacI])
    restriction_analysis = Analysis(rb, transcript_sequence, linear=True)
    data = []
    for oligo_name, insert in zip(transcript_names, inserts):
        ligated_sequence = simulate_ligation_cloning(transcript_sequence, insert, restriction_analysis)
        UTR5, CDS, UTR3 = get_transcript_parts(ligated_sequence)
        data.append({'oligo_name': oligo_name,
            'sequence': ligated_sequence,
            'UTR5': UTR5,
            'CDS': CDS,
            'UTR3': UTR3
        })
    df = pd.DataFrame(data)
    df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    logging.info(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reconstruct viromic MPRA")
    parser.add_argument("--plasmid_gb", required=True, help="Path to the plasmid GenBank file")
    parser.add_argument("--fasta_file", required=True, help="Path to the FASTA file containing oligo sequences")
    parser.add_argument("--output_file", required=True, help="Path to the output TSV file")
    parser.add_argument("--add_sites", action="store_true", help="Flag to add NotI and SacI sites to the inserts")

    args = parser.parse_args()
    main(args.plasmid_gb, args.fasta_file, args.output_file, args.add_sites)