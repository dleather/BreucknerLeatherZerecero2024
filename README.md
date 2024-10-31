# **Bunching in Real Estate Markets: Regulated Building Heights in New York City**

[Link to Paper](https://doi.org/10.1016/j.jue.2024.103683)

[Full Repository w/ Data](https://data.mendeley.com/datasets/3thdvcygww/1)

Blog Post:

## **Overview**

This repository contains the replication materials for the paper:

**Title**: *Bunching in Real Estate Markets: Regulated Building Heights in New York City,* (Journal of Urban Economics, 2024)

**Authors**: Jan K. Brueckner, David Leather, Miguel Zerecero

**Abstract**:

This paper presents a real-estate application of the bunching methodology widely used in other areas of applied microeconomics. The focus is on regulated building heights in New York City, where developers can exceed a parcel’s regulated height by incurring additional costs. Using the bunching methodology, we estimate the magnitude of these extra costs, with the results showing a modest increase in the marginal cost of floor space beyond the regulated building height. We use these estimates to predict the additional floor space that would be created by complete removal of building-height regulation in NYC. While this last exercise is circumscribed by our focus on a limited number of zoning categories, the results suggest that New York could secure notably more housing through lighter height regulation.

---

## **Repository Structure**

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

### **Data Acquisition**

Due to size constraints and licensing issues, raw and processed data files are **not included** in this repository, but can be found in the full repository link above.

### **Instructions**

See `README.pdf`.

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

- **Original Authors**: Jan K. Brueckner, David Leather, Miguel Zerecero
- **Data Sources**:
  - NYC Planning Open Data Portal

---

## **Contact**

For any questions or issues regarding this repository, please contact:

**David Leather**

- **Email**: [david.a.leather@gmail.com](mailto:david.a.leather@gmail.com)
- **GitHub**: [dleather](https://github.com/dleather)

---

### **Notes**

- **Julia Module Details**:
  - The `Bunching.jl` module is essential for running the bunching analysis.
  - It contains custom functions that implement the methodologies discussed in the paper.
  - Ensure that Julia is installed and properly configured to run this module.
