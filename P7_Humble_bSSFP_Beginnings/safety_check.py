import sys
import os
import pypulseq as pp
import matplotlib.pyplot as plt

# Predefined acoustic resonance frequencies (example values - adjust as needed)
DEFAULT_ACOUSTIC_RESONANCES_3T = [ # For our Siemens 3T Prisma
    {"frequency": 590, "bandwidth": 100},
    {"frequency": 1140, "bandwidth": 220}
]

DEFAULT_ACOUSTIC_RESONANCES_7T = [ # For our Siemens 3T Prisma
    {"frequency": 550, "bandwidth": 100},
    {"frequency": 1100, "bandwidth": 300}
]

def extract_md5_hash(seq_file):
    """
    Extract the MD5 hash from the last line of the .seq file
    """
    try:
        with open(seq_file, 'r') as f:
            lines = f.readlines()
            # Look for the hash in the last few lines
            for line in reversed(lines):
                line = line.strip()
                if line.startswith("Hash "):
                    # Extract the hash part after "Hash "
                    hash_value = line[5:]  # Remove "Hash " prefix
                    return hash_value
        return "Hash not found"
    except Exception as e:
        print(f"Error reading seq file for hash: {e}")
        return "Hash extraction failed"

def check_PNS(seq, gradFile, seq_file):
    """
    Plots the PNS level of the sequence using the gradient information in gradFile and throws an error if PNS level is exceeded
    """
    # Check if gradient file exists
    gradient_file = os.path.join(os.path.dirname(__file__), gradFile)
    if os.path.exists(gradient_file):
        print(f"Gradient file found: {gradFile}, using it for PNS check.")
    else:
        print(f"Warning: Gradient file '{gradFile}' not found. Skipping PNS check.")
        return

    # Extract MD5 hash for the title
    md5_hash = extract_md5_hash(seq_file)

    pns_ok, pns_n, pns_c, tpns = seq.calculate_pns(gradient_file, do_plots=True)

    # Add title with MD5 hash
    plt.suptitle(f"PNS Check - Hash: {md5_hash}", fontsize=14, fontweight='bold')

    # Save the plot as PNG
    output_filename = f"pns_check_{os.path.basename(seq_file).replace('.seq', '')}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"PNS plot saved as: {output_filename}")

    if pns_ok:
        print("PNS check passed successfully: maximum PNS level =", 100*pns_n.max(), "%")
        plt.show()
    else:
        print("PNS check failed: maximum PNS level =", 100*pns_n.max(), "%")
        plt.show()
        raise RuntimeError("PNS check failed!")

def check_acoustic_resonances(seq, gradFile, seq_file, scanner):
    """
    Plots the acoustic resonances of the sequence against the forbidden frequencies specified in gradFile
    """
    gradient_file = os.path.join(os.path.dirname(__file__), gradFile) # Path to the Prisma gradient file

    if os.path.exists(gradient_file):
        # Use gradient file if available
        from pypulseq.utils.siemens.readasc import readasc
        asc, extra = readasc(gradient_file) # Parse the ASC file

        # Extract acoustic resonances
        frequencies = asc["asGPAParameters"][0]["sGCParameters"]["aflAcousticResonanceFrequency"]
        bandwidths  = asc["asGPAParameters"][0]["sGCParameters"]["aflAcousticResonanceBandwidth"]
        #print("Frequencies:", frequencies)
        #print("Bandwidths:", bandwidths)

        # Combine frequencies and bandwidths into a list of dictionaries
        acoustic_resonances = [
            {"frequency": freq, "bandwidth": bw}
            for freq, bw in zip(frequencies.values(), bandwidths.values()) # Use .values() to get the dictionary values
            if freq != 0  # Exclude entries with zero frequency
        ]
        print(f"Using acoustic resonances from gradient file: {gradFile}")
    else:
        # Use default acoustic resonances if gradient file is not available
        if scanner == '3T':
            acoustic_resonances = DEFAULT_ACOUSTIC_RESONANCES_3T
        elif scanner == '7T':
            acoustic_resonances = DEFAULT_ACOUSTIC_RESONANCES_7T
        else:
            raise ValueError("Unknown scanner type. Please specify '3T' or '7T'.")
        print(f"Warning: Gradient file '{gradFile}' not found. Using default acoustic resonances.")

    # Extract MD5 hash for the title
    md5_hash = extract_md5_hash(seq_file)

    seq.calculate_gradient_spectrum(acoustic_resonances=acoustic_resonances, use_derivative=False, frequency_oversampling=10, window_width=(min(0.02,seq.duration()[0])))

    # Add title with MD5 hash
    plt.suptitle(f"Acoustic Resonances - Hash: {md5_hash}", fontsize=14, fontweight='bold')

    # Save the plot as PNG
    output_filename = f"acoustic_resonances_{os.path.basename(seq_file).replace('.seq', '')}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Acoustic resonances plot saved as: {output_filename}")
    plt.show()

def safety_check(seq_file, scanner='3T'):
    """
    Test a pulseq sequence for PNS and Acoustic Resonances
    """
    gradFile = "MP_GPA_K2309_2250V_951A_AS82.asc" # Siemens 3T Prisma gradient file for PNS and acoustic resonance check (place in the same folder as this file)

    # Create a Sequence object
    seq = pp.Sequence()

    # Print the sequnece duration
    #print(f"Sequence duration: {seq.duration()[0]:.3f} s")

    # Load the sequence from a .seq file
    seq.read(seq_file)

    # Check the PNS level
    check_PNS(seq, gradFile, seq_file)

    # Check the acoustic resonances
    check_acoustic_resonances(seq, gradFile, seq_file, scanner)

if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("Usage: python safety_check.py <path_to_seq_file> <scanner_type (3T or 7T, default 3T)>")
        sys.exit(1)

    seq_file = sys.argv[1]
    safety_check(seq_file)