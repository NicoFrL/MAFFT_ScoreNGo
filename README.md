# MAFFT_ScoreNGo

A Python tool that screens MAFFT parameters, performs and evaluates alignments, and identifies the optimal alignment strategy given your amino acid sequence dataset.

Created by Nicolas-Frédéric Lipp, PhD. 

## Features
MAFFT_ScoreNGo is useful for determining the theoretically best alignment method using MAFFT** for your specific protein sequence dataset. All you need is an input .fasta file of unaligned sequences. It is recommended to run it on a small dataset - you can use the script I provided, "Rand_NSamp_MyFasta.py" _(available soon)_, to extract 200 sequences or less if your dataset is too large. 

- Automated screening of multiple MAFFT parameters
- Evaluation of alignment quality using custom scoring metrics
- Identification of optimal alignment strategy for given dataset
- Generation of comprehensive result summaries



***MAFFT stands for Multiple sequence Alignment using Fast Fourier Transform. More documentation can be found at [mafft.cbrc.jp](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html), [Katoh et al. (Nucleic Acids Res., 2002)](https://doi.org/10.1093%2Fnar%2Fgkf436), and [Katoh et al. (Brief. Bioinform., 2017)](https://doi.org/10.1093/bib/bbx108).*

## Parameters Overview

MAFFT_ScoreNGo tests various combinations of the following MAFFT parameters:
- Alignment strategies (--genafpair, --localpair, --globalpair)
- Substitution matrices (BLOSUM62, BLOSUM80)
- Gap opening penalties
- Gap extension penalties
- Large gap penalties

For more details about tested parameters and scoring algorithms, please see [PARAMETERS.md](./PARAMETERS.md)


## Installation

1. Clone this repository:  
   ```git clone https://github.com/yourusername/MAFFT_ScoreNGo.git```  
   ```cd MAFFT_ScoreNGo```

2. Install the required Python packages:  
   ```pip3 install -r requirements.txt```

3. Ensure MAFFT is installed and accessible from your command line (see Dependencies section for more details).

## Usage

### Quick Start
1. Run the script with:  
```python3 MAFFT_ScoreNGo.py```

2. Follow the prompts to select your input FASTA file.
Select the screening level : Light, Standard or Aggressive (type ```1```, ```2``` or ```3```)
3. Choose the screening level:
- Light (type `1`): Quick screening with fewer parameter combinations.
- Standard (type `2`): Balanced screening with a moderate number of combinations.
- Aggressive (type `3`): Thorough screening with many parameter combinations.

4. (Optional) Enter your own personalized parameters if desired.

5. Confirm that you want to run the computation on the given number of combinations.

### Customization

You can add custom MAFFT parameters to be tested alongside the predefined combinations. When prompted, enter your parameters in the MAFFT command-line format. For example: ```--maxiterate 1000 --globalpair --thread 4```


### Interpreting Results

The script will output:
1. A ranking of the top 13 alignments based on the final score.
2. Detailed information for each alignment, including parameters used, execution time, and various scores.
3. The best combination of parameters for your dataset.

Results are saved in the `mafft_results` directory.

## Troubleshooting

- If MAFFT is not found, ensure it's correctly installed and added to your system PATH.
- For memory issues with large datasets, try using a smaller subset of sequences.
- If you encounter Python-related errors, verify that all dependencies are correctly installed.

## Dependencies

- Python 3.x (3.7 or later recommended)
- Biopython 1.81
- Tkinter
- MAFFT (must be installed separately and available in your system PATH)

To check if MAFFT is properly installed and available, run the following command in your terminal:  
```mafft --version```

This should display the installed version of MAFFT. If you see an error message instead, please refer to the MAFFT Installation section below.

### Python and Biopython

Ensure you have Python 3.7 or later installed. You can install Biopython and other Python dependencies using:  
```pip3 install -r requirements.txt```

### Tkinter Installation

Tkinter is usually included with Python, but on some systems, it may need to be installed separately:

- On macOS, if you've installed Python via Homebrew:  
 ```brew install python-tk```

- On Ubuntu/Debian Linux:  
  ```sudo apt-get install python3-tk```

- On other systems, please refer to your system's package manager or Python distribution instructions.

### MAFFT Installation

MAFFT must be installed separately and be available in your system PATH. Installation instructions vary by operating system:

- On macOS with [Homebrew](https://brew.sh):  
  ```brew install mafft```

- On Ubuntu/Debian Linux:  
  ```sudo apt-get install mafft```

- For other systems or for manual installation, please refer to the MAFFT official website: https://mafft.cbrc.jp/alignment/software/

After installation, verify MAFFT is accessible by running:  
```mafft --version```

## Performances
MAFFT_ScoreNGo.py was tested on a dataset of 200 sequences with an average length of 183 residues, using MacOS Sonoma on an ARM M1 Max (10-core CPU, 32-core GPU) with 32 GB of RAM. With this configuration, it took around 30 minutes to run all alignments and scoring assessments.
Don't forget to "[Caffeinate](https://www.theapplegeek.co.uk/blog/caffeinate)" your Mac! (or [systemd-inhibit](https://evanhahn.com/systemd-inhibit-alternative-to-macos-caffeinate/) your Linux machine). 



## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have any questions, please open an issue on the GitHub repository.

## Author

Nicolas-Frédéric Lipp, PhD  
https://github.com/NicoFrL

## Acknowledgements

This project was developed with the assistance of AI language models, which provided guidance on code structure, best practices, and documentation. The core algorithm and scientific approach were designed and implemented by the author.
