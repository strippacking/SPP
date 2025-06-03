This README file was generated on 2025-06-03 by Duong Quy Le

======================

GENERAL INFORMATION

*Dataset Title: Benchmark Dataset for Two-Dimensional Strip Packing Problem: 41 Standard Instances from Established Test Sets

Dataset DOI: <TO BE ASSIGNED BY SCIELO DATA>

*Principal Investigator/Author Information
Name: Tuyen Kieu Van
ORCID: 
Institution: VNU University of Engineering and Technology, Hanoi, Vietnam
Email: tuyenkv@vnu.edu.vn

Co-Investigator/Author Information
Name: Duong Quy Le
ORCID: 
Institution: VNU University of Engineering and Technology, Hanoi, Vietnam
Email: 21020560@vnu.edu.vn

Alternate Contact Information
Name: Khanh Van To
ORCID: 
Institution: VNU University of Engineering and Technology, Hanoi, Vietnam
Email: khanhtv@vnu.edu.vn

*Data collection date: 2024-01-01 to 2024-12-31

Geographic location of data collection: Hanoi, Vietnam

Funding sources supporting data collection: VNU University of Engineering and Technology under project number CN24.10

Keywords: Two-dimensional strip packing, benchmark instances, combinatorial optimization, SAT solving, MaxSAT, rectangle packing, cutting and packing problems

===========================

DATA & FILE OVERVIEW

*List of Files:
- ins-1.txt to ins-41.txt: 41 benchmark instances for the Two-Dimensional Strip Packing Problem
- check.py: Python script for validating instance format and basic statistics

Relationship between files, if important:
All instance files (ins-1.txt to ins-41.txt) follow the same format and represent different problem instances from established benchmark datasets in the strip packing literature. The instances are numbered sequentially but correspond to different original dataset sources (HT, BENG, CGCUT, GCUT, NGCUT).

File sizes and formats:
- Instance files: Plain text format (.txt), sizes range from 0.1 KB to 2 KB
- Total dataset size: Approximately 50 KB
- Format: ASCII text with integer values

Related data not included in this data package:
- Optimal solution visualizations (available in the main repository)
- Solver implementation code (available in the main repository)
- Experimental results (SPP_Result.xlsx in main repository)

===========================

METHODOLOGICAL INFORMATION

*Description of methods used for data collection/generation:
The benchmark instances were collected from established datasets in the Two-Dimensional Strip Packing Problem literature. These instances have been widely used in academic research and represent standard test cases for algorithm evaluation. Sources include:
1. HT instances (Hopper & Turton, 2001)
2. BENG instances (Bengtsson, 1982)
3. CGCUT instances (Christofides & Whitlock, 1977)
4. GCUT instances (Beasley, 1985)
5. NGCUT instances (Beasley, 1985)

*Methods for data processing:
The original instances were converted to a standardized format where:
- First line contains the strip width (W)
- Second line contains the number of rectangles (n)
- Following lines contain width and height of each rectangle (wi hi)
All instances were verified for format consistency and correctness.

*Specific instrument or software information needed to interpret the data:
- Any text editor for viewing instance files
- Python 3.x for running the validation script (check.py)
- No specialized software required for basic data interpretation

Standards and calibration information, if appropriate:
All instances follow the standard format for 2D Strip Packing Problem representation as commonly used in the Operations Research community.

Environmental/experimental conditions:
Not applicable - these are mathematical problem instances.

Describe any quality assurance procedures performed on the data:
1. Format validation: All instances checked for correct format compliance
2. Data integrity: Verified that all rectangle dimensions are positive integers
3. Consistency check: Ensured strip widths are appropriate for the given rectangles
4. Literature verification: Cross-referenced with original sources where possible

People involved with sample collection, processing, analysis, and/or submission:
- Tuyen Kieu Van: Dataset compilation and validation
- Duong Quy Le: Format standardization and quality assurance
- Khanh Van To: Literature review and source verification

===========================

DATA-SPECIFIC INFORMATION FOR: Instance Files (ins-1.txt to ins-41.txt)

Variables:
- W (Line 1): Strip width - Integer value representing the fixed width of the strip - Units: abstract units
- n (Line 2): Number of rectangles - Integer value (range: 7 to 200 in this dataset) - Units: count
- wi (Lines 3 to n+2, first value): Width of rectangle i - Positive integer - Units: abstract units  
- hi (Lines 3 to n+2, second value): Height of rectangle i - Positive integer - Units: abstract units

Missing data codes: None - all data is complete

Specialised formats or other abbreviations used:
- Instance naming: ins-X.txt where X ranges from 1 to 41
- Original instance names are mapped as follows:
  * ins-1 to ins-12: HT01 to HT12 (Hopper-Turton instances)
  * ins-13 to ins-15: CGCUT01 to CGCUT03
  * ins-16 to ins-19: GCUT01 to GCUT04  
  * ins-20 to ins-31: NGCUT01 to NGCUT12
  * ins-32 to ins-41: BENG01 to BENG10

DATA-SPECIFIC INFORMATION FOR: check.py

Description: Python validation script that:
- Reads all instance files
- Validates format correctness
- Computes basic statistics (number of rectangles, strip width, etc.)
- Generates summary information about the dataset

Requirements: Python 3.x with standard libraries

==========================

SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: CC BY 4.0 (Creative Commons Attribution 4.0 International)

Access conditions: Open access - no restrictions. Data can be freely used for research and educational purposes with proper attribution.

Links to publications that cite or use the data:
- Van, T.K., Le, D.Q., To, K.V. "Efficient SAT and MaxSAT techniques for solving the Two-Dimensional Strip Packing Problem" (submitted to Pesquisa Operacional, 2024)

Links to other publicly accessible locations of the data:
- GitHub repository: https://github.com/strippacking/SPP (if applicable)

Links/relationships to auxiliary datasets:
- Complete experimental results: SPP_Result.xlsx (available in main repository)
- Solver implementations: Various Python scripts in main repository

Were the data derived from another source?
Yes, the instances were derived from multiple established sources:
1. University of Bologna OR Group library
2. Original research papers cited in the methodological section
3. Standard benchmark collections in the Operations Research community

Recommended citation for this dataset:
Van, Tuyen Kieu; Le, Duong Quy; To, Khanh Van. 2024. Replication data for: Efficient SAT and MaxSAT techniques for solving the Two-Dimensional Strip Packing Problem. SciELO Data, V1.

===========================

DOCUMENTATION & SUPPORT

Links to additional documentation:
- Main repository README.md: Contains comprehensive documentation of the research project
- Research paper: Detailed methodology and experimental analysis
- Instance format specification: Standard 2D Strip Packing Problem format

Support information:
For questions about the dataset or its usage, contact:
- Duong Quy Le: 21020560@vnu.edu.vn
- Tuyen Kieu Van: tuyenkv@vnu.edu.vn
Technical support for reproducing experiments is available through the main repository.

===========================

PRIVACY & ETHICS INFORMATION

Ethical approvals obtained for data collection: Not applicable - this dataset consists of mathematical problem instances from published literature, not human subject data.

Privacy considerations: Not applicable - no personal or sensitive data is included in this dataset.

===========================

PARTICIPANT CONSENT

Details about informed consent: Not applicable - this research does not involve human participants. All data consists of mathematical problem instances from established academic sources.
