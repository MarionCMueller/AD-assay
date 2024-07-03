import sys
import subprocess
import os

def call_vcf(reference_assembly, bamfile_parent1, bamfile_parent2, output_vcf):
    # Check if the output file already exists
    if os.path.exists(output_vcf):
        print(f"Output file {output_vcf} already exists. Skipping VCF calling for this sample.")
        return

    call_command = f"freebayes -f {reference_assembly} {bamfile_parent1} {bamfile_parent2} --haplotype-length 0 --min-alternate-count 20 --min-alternate-fraction 0 --pooled-continuous --limit-coverage 400 > {output_vcf}"
    
    print("Executing command:", call_command)
    
    try:
        subprocess.run(call_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error:", e)
        sys.exit(1)

def parse_vcf_positions(input_vcf):
    positions = []
    with open(input_vcf, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 5 and fields[3] != "." and fields[4] != ".":
                positions.append((fields[0], int(fields[1])))
    return positions


def filter_vcf(input_vcf, output_txt, min_percentage=95, min_dp=10):
    with open(input_vcf, 'r') as vcf_file, open(output_txt, 'w') as output_file:
        for line_number, line in enumerate(vcf_file, start=1):
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # Ensure there are enough elements in the line
            if len(fields) < 10:
                print(f"Error: Line {line_number} does not have enough elements: {line.strip()}. Skipping this line.")
                continue

            if fields[9] == "." or fields[10] == ".":
                continue

            # Skip non-SNPs or positions with more than two alleles
            if len(fields[4].split(",")) == 1: 
               

              chrom = fields[0]
              pos = int(fields[1])

             # Get genotype information for the sample
              genotype_info_parent1 = fields[9].split(":")
              genotype_info_parent2 = fields[10].split(":")

              allelecount_parent1 = genotype_info_parent1[2]
              allelecount_parent2 = genotype_info_parent2[2]

              if "," in allelecount_parent1:

               allelecount_parent_split1 = allelecount_parent1.split(",")
               total_reads_parent1 = int(allelecount_parent_split1[0]) + int(allelecount_parent_split1[1])
               if total_reads_parent1 != 0:
                ref_al_prop1 = int(allelecount_parent_split1[0])/total_reads_parent1
                alt_al_prop1 = int(allelecount_parent_split1[1])/total_reads_parent1
              elif "," not in allelecount_parent1:
               allelecount_parent_split1 = allelecount_parent1
               total_reads_parent1 = int(allelecount_parent_split1) 
               ref_al_prop1 = 1
               alt_al_prop1 = 0
              if "," in allelecount_parent2:

               allelecount_parent_split2 = allelecount_parent2.split(",")
               total_reads_parent2 = int(allelecount_parent_split2[0]) + int(allelecount_parent_split2[1])
               if total_reads_parent2 != 0:
                ref_al_prop2 = int(allelecount_parent_split2[0])/total_reads_parent2
                alt_al_prop2 = int(allelecount_parent_split2[1])/total_reads_parent2
              elif "," not in allelecount_parent2:
               allelecount_parent_split2 = allelecount_parent2
               total_reads_parent2 = int(allelecount_parent_split2) 
               ref_al_prop2 = 1
               alt_al_prop2 = 0
             
            


             
            
            if genotype_info_parent1[0] == '0/0' and genotype_info_parent2[0] == '1/1' and ref_al_prop1 >= 0.95 and alt_al_prop2 >= 0.95 and total_reads_parent1 > 10 and total_reads_parent2 > 10: 
               pos1 = int(pos)-1
               output_file.write(f"{chrom}\t{pos1}\t{pos}\n")
               
            	# Check if the sample is homozygous for the alternative allele
            elif genotype_info_parent1[0] == '1/1' and genotype_info_parent2[0] == '0/0' and alt_al_prop1 >= 0.95 and ref_al_prop2 >= 0.95 and total_reads_parent1 > 10 and total_reads_parent2 > 10:
               pos1 = int(pos)-1
               output_file.write(f"{chrom}\t{pos1}\t{pos}\n")
               

if len(sys.argv) != 2:
    print("Usage: python script.py input_file.txt")
    sys.exit(1)
    

def call_vcf_filtered(reference_assembly, positions_file, bamfile_pool, output_vcf):
    # Check if the output file already exists
    if os.path.exists(output_vcf):
        print(f"Output file {output_vcf} already exists. Skipping SNP calling for this sample.")
        return

    freebayes_command = f"freebayes -f {reference_assembly}  {bamfile_pool} --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic  -t {positions_file} > {output_vcf} "
    
    try:
        subprocess.run(freebayes_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error:", e)
        sys.exit(1)


def parse_final_snp_positions(input_vcf, output_txt):
    with open(input_vcf, 'r') as vcf_file, open(output_txt, 'w') as output_file:
        # Write header to the output file
        output_file.write("Chromosome\tPosition\tCount_reference_allele\tCount_alternative_allele\n")

        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            try:
                chrom = fields[0]
                pos = int(fields[1])
                format_info = fields[9].split(":")
                if len(fields[4].split(",")) == 1: 

                 # Extract counts from the genotype information
                 counts = format_info[2]
                 print(f"{counts}")
                
                 # Ensure the counts are present and in the expected format
                 if "," in format_info[2]: 
                  counts = format_info[2].split(",")
                  count_ref_allele = int(counts[0]) 
                  count_alt_allele = int(counts[1]) 

                  output_file.write(f"{chrom}\t{pos}\t{count_ref_allele}\t{count_alt_allele}\n")
                 elif "," not in format_info[2]: 
                  count_ref_allele = int(format_info[2]) 
                  count_alt_allele = 0
                  output_file.write(f"{chrom}\t{pos}\t{count_ref_allele}\t{count_alt_allele}\n")

            except IndexError as e:
                print(f"Error: Index error in line {line.strip()}. Skipping this line. Details: {e}")
                continue
            except ValueError as e:
                print(f"Error: Value error in line {line.strip()}. Skipping this line. Details: {e}")
                continue



# Input file with data for each sample
input_file = sys.argv[1]

# Output directory for VCF files
output_vcf_dir = "output.v8"

# Output text file with final filtered SNP positions
output_txt_final = "final_filtered_positions.txt"

# Read data from the input file and perform VCF calling and filtering
with open(input_file, 'r') as input_data:
    next(input_data)  # Skip header line
    for line in input_data:
        fields = line.strip().split("\t")
        bamfile_pool, name_pool, bamfile_parent1, name_parent1, bamfile_parent2, name_parent2, reference_assembly, selection_line = fields

        # Step 1: Create output directory if it doesn't exist and perform bcftools mpileup and bcftools call for parents
        os.makedirs(output_vcf_dir, exist_ok=True)
        output_vcf1_path = os.path.join(output_vcf_dir, f"{name_parent1}.{name_parent2}.{reference_assembly}.vcf")
        output_txt_parent1 = os.path.join(output_vcf_dir, f"{name_parent1}.{name_parent2}.{reference_assembly}.filtered_positions.txt")

        call_vcf(reference_assembly, bamfile_parent1, bamfile_parent2, output_vcf1_path)

        # Check if the output file from step 1 exists before proceeding to step 2
        if not os.path.exists(output_vcf1_path):
            print(f"Error: Output file {output_vcf1_path} does not exist. Skipping to the next iteration.")
            continue

        # Step 2: Apply filtering and write the results to the output file
        filter_vcf(output_vcf1_path, output_txt_parent1, min_percentage=95, min_dp=10)

        # Perform bcftools mpileup and bcftools call using filtered positions for the pool
        output_vcf2_path = os.path.join(output_vcf_dir, f"{name_pool}.{reference_assembly}.{selection_line}.vcf")
        output_txt_filtered_path = os.path.join(output_vcf_dir, f"{name_parent1}.{name_parent2}.{reference_assembly}.filtered_positions.txt")
        
        # Call bcftools mpileup and bcftools call using filtered positions for the pool
        call_vcf_filtered(reference_assembly, output_txt_filtered_path, bamfile_pool, output_vcf2_path)
        print(f"Filtering and additional SNP calling completed successfully for {name_parent1}_{name_parent2} and {name_pool}.")
        
        # Step 4: Parse the final SNP positions and create the output file
        output_final_snp_positions = os.path.join(output_vcf_dir, f"{name_pool}.{reference_assembly}.{selection_line}.final_snp_positions.txt")
        parse_final_snp_positions(output_vcf2_path, output_final_snp_positions)



print("Script execution completed.")

