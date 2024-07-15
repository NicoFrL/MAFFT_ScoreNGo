# MAFFT_ScoreNGo

A Python tool that screens MAFFT parameters, performs and evaluates alignments, and identifies the optimal alignment strategy given your amino acid sequence dataset.

Created by Nicolas-Frédéric Lipp, PhD. 

## Installation

1. Clone this repository:
   git clone https://github.com/yourusername/MAFFT_ScoreNGo.git
   cd MAFFT_ScoreNGo

2. Install the required Python packages:
   pip3 install -r requirements.txt

3. Ensure MAFFT is installed and accessible from your command line (see Dependencies section for more details).

## Usage

Run the script with:
python3 MAFFT_ScoreNGo.py

Follow the prompts to select your input FASTA file.

## Output

The script will create a directory named 'mafft_results' in the same location as your input file, containing:
- Individual alignment files
- A summary of results (mafft_results_summary.txt)
- A list of MAFFT commands used (mafft_commands.txt)

## Dependencies

- Python 3.x (3.7 or later recommended)
- Biopython 1.81
- Tkinter
- MAFFT (must be installed separately and available in your system PATH)

To check if MAFFT is properly installed and available, run the following command in your terminal:
mafft --version

This should display the installed version of MAFFT. If you see an error message instead, please refer to the MAFFT Installation section below.

### Python and Biopython

Ensure you have Python 3.7 or later installed. You can install Biopython and other Python dependencies using:
pip3 install -r requirements.txt

### Tkinter Installation

Tkinter is usually included with Python, but on some systems, it may need to be installed separately:

- On macOS, if you've installed Python via Homebrew:
  brew install python-tk

- On Ubuntu/Debian Linux:
  sudo apt-get install python3-tk

- On other systems, please refer to your system's package manager or Python distribution instructions.

### MAFFT Installation

MAFFT must be installed separately and be available in your system PATH. Installation instructions vary by operating system:

- On macOS with [Homebrew](https://brew.sh):  
  brew install mafft

- On Ubuntu/Debian Linux:  
  sudo apt-get install mafft

- For other systems or for manual installation, please refer to the MAFFT official website: https://mafft.cbrc.jp/alignment/software/

After installation, verify MAFFT is accessible by running:
mafft --version

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
