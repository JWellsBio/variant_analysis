# variant_analysis
## The purpose of the project is to determine how effective liquid biopsies can be in detecting metastases by detecting variants found in the DNA in plasma

## Samples
### Samples include:
#### Patient 10:
Liver 1 (met)
Liver 2A (met)
Liver 5 (met)
Plasma

#### Patient 9:
Lymph Node (met)
Omental (met)
Ovary (met)
Plasma

#### Patient 2
2 Liver mets
2 Breast primaries
*Normal Liver instead of buffy*

#### Patient 8
Axillary (met)
2 Breast primaries
Plasma

#### Patient 9A
Heart (met)
Left Kidney (met)
Right Kidney (met)
2 Liver mets
2 Omental mets
Plasma
*Uses patient 9 buffy as normal*

### Filters used on variant calls
- Indels removed
- Intron and intergenic variants removed
- All Phred-scale quality measures > 20 (99.9% accuracy)
- Plasma calls must have MAF > 0.01
- Tumor calls must have MAF > 0.05
- All Mutect2 calls must also be similarly called by Varscan or Freebayes

