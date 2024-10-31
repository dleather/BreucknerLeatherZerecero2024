ing in Real Estate Markets: Regulated Building Heights in New York City**

## **Overview**

This repository contains the replication materials for the paper:

**Title**: *Bunching in Real Estate Markets: Regulated Building Heights in New York City*

**Authors**: Jan K. Brueckner, **David Leather**, Miguel Zerecero

**Abstract**:

This paper presents a real-estate application of the bunching methodology widely used in other areas of applied microeconomics. The focus is on regulated building heights in New York City, where developers can exceed a parcel’s regulated height by incurring additional costs. Using the bunching methodology, we estimate the magnitude of these extra costs, with the results showing a modest increase in the marginal cost of floor space beyond the regulated building height. We use these estimates to predict the additional floor space that would be created by complete removal of building-height regulation in NYC. While this last exercise is circumscribed by our focus on a limited number of zoning categories, the results suggest that New York could secure notably more housing through lighter height regulation.

---

## **Repository Structure**

itional_stata_code/ ├── Bunching/ │ └── src/ │ └── Bunching.jl ├── processed_data/ ├── raw_data/ │ ├── MAPPLUTO/ │ └── PLUTO/ │ └── extra_csv_files/ ├── results/ │ ├── figures/ │ └── tables/ └── scripts/

ditional_stata_code/`**: Additional Stata scripts used in the analysis.
- **`scripts/`**: Main scripts for data processing and analysis.
- **`Bunching/src/Bunching.jl`**: A Julia module containing all custom functions used in the bunching methodology.
- **`raw_data/`**: Contains raw datasets (**not included in this repository**).
  - **`MAPPLUTO/`**
  - **`PLUTO/`**
- **`processed_data/`**: Processed datasets (**not included in this repository**).
- **`results/`**:
  - **`figures/`**: Output figures from the analysis.
  - **`tables/`**: Output tables from the analysis.

---

## **Getting Started**

### **Prerequisites**

- **Software Requirements**:
  - Stata 15 or higher
  - Julia (version compatible with the code)
  - Quarto (for rendering documents)
  - Other software as required (e.g., Python, R)

### **Data Acquisition**

Due to size constraints and licensing issues, raw and processed data files are **not included** in this repository. Please follow the instructions below to obtain the necessary data:

1. **MAPPLUTO Data**:
   - Download from the [NYC Planning Open Data Portal](https://www1.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page).
   - Place the downloaded files in `raw_data/MAPPLUTO/`.

2. **PLUTO Data**:
   - Also available from the [NYC Planning Open Data Portal](https://www1.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page).
   - Place the files in `raw_data/PLUTO/`.

### **Running the Analysis**

1. **Data Processing**:
   - Navigate to the `scripts/` directory.
   - Run the data cleaning and processing scripts.

     ```stata
     do data_processing.do
     ```

2. **Bunching Analysis**:
   - Navigate to `Bunching/src/` and run the Julia module `Bunching.jl`.

     ```julia
     # In Julia REPL or script
     include("Bunching/src/Bunching.jl")
     ```

   - The `Bunching.jl` module contains all custom functions used in the bunching methodology analysis.

3. **Generating Results**:
   - The analysis scripts will generate figures and tables.
   - Outputs are saved in the `results/` directory.

### **Reproducing Figures and Tables**

- Ensure all scripts are executed in the correct order.
- Update any file paths in the scripts if necessary.
- Refer to the comments within the scripts and the Julia module for additional guidance.

---

## **Project Description**

This project replicates the analysis conducted in the paper "*Bunching in Real Estate Markets: Regulated Building Heights in New York City*". The study applies the bunching methodology to examine how regulated building heights impact the real estate market in NYC.

### **Key Objectives**

- **Estimate Extra Costs**: Quantify the additional costs developers incur when exceeding regulated building heights.
- **Predict Additional Floor Space**: Assess the potential increase in floor space from removing building-height regulations.
- **Policy Implications**: Provide insights into how lighter height regulations could alleviate housing shortages in NYC.

---

## **Results**

The replication study confirms the findings of the original paper:

- **Marginal Cost Increase**: There is a modest increase in the marginal cost of floor space beyond the regulated building height.
- **Potential Housing Increase**: Easing height regulations could lead to a notable increase in available housing in New York City.
- **Economic Implications**: The results suggest that policy adjustments in building regulations can have significant effects on urban development.

---

## **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## **Acknowledgments**

- **Original Authors**: Jan K. Brueckner, **David Leather**, Miguel Zerecero
- **Data Sources**:
  - NYC Planning Open Data Portal
- **Contributors**:
  - Collaborators and reviewers who assisted in the replication study.

---

## **Contact**

For any questions or issues regarding this repository, please contact:

**David Leather**

- **Email**: [your.email@example.com](mailto:your.email@example.com)
- **GitHub**: [your-github-username](https://github.com/your-github-username)

---

## **Additional Resources**

- **Original Paper**: [Link to the published paper if available online]
- **Quarto Blog Post**: [Link to the blog post if you've created one]
- **Related Works**: References to other studies or papers on related topics.

---

### **Notes**

- **Julia Module Details**:
  - The `Bunching.jl` module is essential for running the bunching analysis.
  - It contains custom functions that implement the methodologies discussed in the paper.
  - Ensure that Julia is installed and properly configured to run this module.
