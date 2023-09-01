import pysam

# Define input and output BAM file paths
input_bam_path = "HG00171.mapped.ILLUMINA.bwa.FIN.exome.20120522.bam" # path to file
output_bam_path = "output.bam"

# Define SNV information
chromosome = "1"
position = 77876703  # The position where you want to introduce the SNV (1-based)
new_base = "C"   # The new base to replace the reference base

# Open input BAM file for reading
input_bam = pysam.AlignmentFile(input_bam_path, "rb")

# Create an output BAM file for writing
output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=input_bam.header)

# Iterate through the reads in the input BAM file: select first chromosome and small region
for read in input_bam.fetch('1', 77876700,77876705): 
    try: 
        print(read.reference_start)
        # Check if the read overlaps with the SNV position
        if read.reference_name == chromosome and read.reference_start <= position <= read.reference_end:
            # Modify the sequence and quality scores to introduce the SNV
            read_sequence = list(read.seq)
            read_quality = list(read.qual)
    
            # Calculate the position within the read
            read_pos = position - read.reference_start-1
            print(read_pos, len(read.seq), read.seq)
    
            # Replace the reference base with the new base
            read_sequence[read_pos] = new_base
            read_quality[read_pos] = chr(ord(read_quality[read_pos]) + 1)  # Increment quality score
    
            # Update the read's sequence and quality scores
            read.seq = "".join(read_sequence)
            read.qual = "".join(read_quality)
    except TypeError:
        continue
    # Write the modified read to the output BAM file
    output_bam.write(read)

# Close the BAM files
input_bam.close()
output_bam.close()
