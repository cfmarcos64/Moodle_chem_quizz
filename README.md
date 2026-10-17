# 🧪 Moodle Chemistry Suite

![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-green)
![Streamlit](https://img.shields.io/badge/Streamlit-Cloud-orange)

A comprehensive **Streamlit** suite designed to assist chemistry educators in quickly, intuitively, and professionally generating Moodle-compatible XML quiz banks (short answer and multiple choice questions with embedded molecular structure images).

All tools share a unified design philosophy: the educator provides compound names or SMILES strings, and the suite handles automatic NCI CIR name-to-SMILES resolution, stereochemical analysis via **RDKit**, high-quality molecular rendering, and direct XML export ready for importing straight into Moodle.

---

## 🚀 Included Modules

1. **Systematic Nomenclature (`Nomenclature_for_Moodle.py`)**  
   * Generates short answer questions evaluating IUPAC systematic nomenclature for organic compounds.
   * Features automated validation and support for multiple accepted synonym answers.

2. **Functional Groups (`FG_for_Moodle.py`)**  
   * Identifies and analyses functional groups present within molecular structures.
   * Generates visually highlighted structures for short answer identification quizzes.

3. **Stereochemistry Determination (`Stereochemistry_for_Moodle.py`)**  
   * Automatically detects chiral centres ($R/S$) and geometric double bonds ($Z/E$) using RDKit.
   * Dynamically adapts question prompts according to the detected stereogenic element and exports to Moodle short answer format.

4. **Stereoisomer Comparison & Analysis (`Compare_for_Moodle.py`)**  
   * Evaluates the stereochemical relationship between a base molecule and either automatically enumerated stereoisomers or a second manually entered compound.
   * Rigorously classifies pairs into four distinct categories: Enantiomers, Diastereomers, Identical compounds (including meso forms), and Different non-stereoisomeric compounds (exporting as Moodle multiple choice questions).

---

## ✨ Key Features

- **Bilingual Interface**: Full support for both English and Spanish with instant language switching in the sidebar.
- **Robust Workflows**: Supports individual entry (with full Return/Enter key submission) as well as bulk processing via Excel (`.xlsx`) or CSV files.
- **Standardised UI/UX**: Professional two-column layout, clear tabbed forms, persistent thumbnail previews, and a clean, responsive workflow.
- **Advanced Chemical Intelligence**: Powered by RDKit for robust stereochemical calculations and NCI CIR for seamless name-to-SMILES resolution.
- **Optimised XML Export**: Prominent green export buttons positioned at the top of the question queue for quick download and effortless Moodle integration.

---

## 🛠️ How to Run the Suite

There are two primary ways to access and use this tool:

### Option 1: Use the Public Web Application

The application is deployed publicly and can be accessed directly through this link:

👉 https://chemformoodle.streamlit.app/

### Option 2: Run Locally 

1. **Clone the Repository & Install Python Dependencies**:

```bash
git clone https://github.com/cfmarcos64/Moodle_chem_quizz.git
cd Moodle_chem_quizz

pip install -r requirements.txt
```

## 📦 Project Dependencies

### `requirements.txt`
```text
streamlit>=1.31.0
rdkit>=2023.09.1
pandas>=2.0
openpyxl>=3.1
requests>=2.28
pillow>=10.0
```

### `packages.txt` (For Linux / Streamlit Cloud RDKit Stability)
```text
libxrender1
libxext6
```

---

## 📁 Repository Structure

The repository is organised into a clean, modular structure:

```text
Moodle_chem_quizz/
├── Chem_quiz.py                   # Master navigation script & global language state
├── Nomenclature_for_Moodle.py     # Systematic Nomenclature module
├── FG_for_Moodle.py               # Functional Groups module
├── Stereochemistry_for_Moodle.py  # Stereochemistry module (R/S and Z/E)
├── Compare_for_Moodle.py          # Stereoisomer Comparison module
├── requirements.txt               # Python package dependencies
├── packages.txt                   # Linux system dependencies for RDKit
└── README.md                      # Project documentation
```

---

## 📄 Licence

This project is licensed under the CC BY-NC-SA 4.0 Licence. See the [LICENSE](LICENSE) file for further details.

---

## 👤 Author & Contact
* **Created by**: Carlos Fernández Marcos
* **Institution**: Universidad de Extremadura (University of Extremadura, Spain)
* **GitHub**: [@cfmarcos64](https://github.com/cfmarcos64)
* **Email**: cfmarcos@unex.es
